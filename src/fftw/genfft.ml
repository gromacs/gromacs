(*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *  
 *)

(* $Id$ *)

#open "string";;
#open "printf";;

let CVSId = "$Id$";;

(* the expressions our program manipulates *)
type Expr = Var of string 
  | Subscript of Expr * Expr   (* array subscripts *)
  | Real of float 
  | Integer of int 
  | Plus of Expr list 
  | Times of Expr * Expr 
  | Uminus of Expr
  | Comma of Expr * Expr
  | Assign of Expr * Expr
  | Call of Expr * Expr
  | Binop of string * Expr * Expr
;;

let real = "FFTW_REAL" ;;
let precision_format = "%.16g";;

(* conversion from float to strings *)
let string_of_float f =
  let s = format_float precision_format f in
  try
    for i = 0 to pred(string_length s) do
      match nth_char s i with `.` | `e` | `E` -> raise Exit | _ -> ()
    done;
    s ^ ".0"
  with Exit ->
    s
;;

(* conversion to float to konst *)
let konst_of_float f =
  let f2 = if (f >=. 1.0) then (f -. (float_of_int (int_of_float f))) else f
  in let q = string_of_int (int_of_float(f2 *. 1.0E9))
  in let r = "0000000000" ^ q
  in let l = string_length r
  in
      if (f >=. 1.0) then
	("FFTW_K" ^ (string_of_int (int_of_float f)) ^ "_" ^ 
         (sub_string r (l - 9) 9))
      else
      	("FFTW_K" ^ (sub_string r (l - 9) 9))
;;

(* define this to be either konst_of_float or string_of_float *)
let c_code_of_float = konst_of_float;;

(* auxiliary functions *)
let rec foldr f a = function
      [] -> a
    | x::y -> f x (foldr f a y)
;;
    
let rec forall combiner a b f =
    if (a >= b) then []
    else combiner (f a) (forall combiner (a + 1) b f)
;;

let all pred l = foldr (prefix &&) true (map pred l);;
let exists pred l = foldr (prefix or) false (map pred l);;
(* let flatten = foldr (prefix @) [];; *)
let compose f g = fun x -> f (g x);;

(* Hmm... ML won't allow second-order polymorphism.  Oh well.. *)
(* let forall_flat = forall (prefix @);; *)
let rec forall_flat a b f = 
    if (a >= b) then []
    else (f a) @ (forall_flat (a + 1) b f)
;;

(* 
 * the unparser transforms expressions into strings 
 *)
let rec unparse_expr =
    let rec unparse_plus = function
	     [] -> ""
	     | (Uminus a :: b) -> " - " ^ (parenthesize a) ^ (unparse_plus b)
	     | (a :: b) -> " + " ^ (parenthesize a) ^ (unparse_plus b)
    and parenthesize x = match x with
	(Integer i) -> unparse_expr x
      |	(Real f) -> unparse_expr x
      |	(Var v) -> unparse_expr x
      |	(Subscript (a,b)) -> unparse_expr x
      |	_ -> "(" ^ (unparse_expr x) ^ ")"

    in function
	 Var x -> x
       | Subscript (a, i) -> (unparse_expr a) ^ "[" ^
	                     (unparse_expr i) ^ "]"
       | Real f -> "((" ^ real ^ ") " ^ c_code_of_float f ^ ")"
       | Integer i -> string_of_int i
       | Plus [] -> "0.0 /* error */"
       | Plus [a] -> (unparse_expr a)
       | Plus (a::b) -> (parenthesize a) ^ (unparse_plus b)
       | Times (a, b) -> (parenthesize a) ^ " * " ^ (parenthesize b)
       | Comma (a, b) -> (unparse_expr a) ^ ", " ^ (unparse_expr b)
       | Assign (a, b) -> (unparse_expr a) ^ " = " ^ (unparse_expr b)
       | Call (a, b) -> (unparse_expr a) ^ "(" ^ (unparse_expr b) ^ ")"
       | Uminus a -> "(- (" ^ (unparse_expr a) ^ "))"
       | Binop (op, a, b) -> (unparse_expr a) ^ op ^ (unparse_expr b)
;;


(********************************************
      Expression simplifier
 ********************************************)
let epsilon = 1.0E-12;;

(* comparison predicate for real numbers *)
let almost_equal x y = (abs_float (x -. y) <. epsilon);;

(* some expressions are simple, some are not *)
let rec simple = function
	 Var c -> true
       | Real f -> true
       | Integer f -> true
       | Plus a -> false
       | Times (a, b) -> false
       | Uminus e -> simple e
       | Subscript (a, b) -> false
       | Call (a, b) -> simple b
       | _ -> false
;;

let rec simplify = 
  let rec reduce_sum x = match x with
	 [] -> []
       | (Plus a) :: b -> reduce_sum (a @ b)    (* flatten *)
       | (Uminus (Plus a)) :: b -> 
             reduce_sum ((map (fun x -> (simplify (Uminus x))) a) @ b)
       | [Real a] -> if (almost_equal a 0.0) then [] else x
       | [Uminus (Real a)] -> if (almost_equal a 0.0) then [] else x
       | [Integer a] -> if (a == 0) then [] else [Integer a]
       | (Real a) :: (Real b) :: s -> reduce_sum (Real (a +. b) :: s)
       | (Uminus (Real a)) :: (Real b) :: s -> 
					 reduce_sum (Real (-. a +. b) :: s)
       | (Real a) :: (Uminus (Real b)) :: s -> 
					 reduce_sum (Real (a -. b) :: s)
       | (Uminus (Real a)) :: (Uminus (Real b)) :: s -> 
					 reduce_sum (Real (-. a -. b) :: s)
       | (Integer a) :: (Integer b) :: s ->
				      reduce_sum (Integer (a + b) :: s)
       | (Real a) :: b :: s -> b :: reduce_sum ((Real a) :: s)
       | (Uminus (Real a)) :: b :: s -> b :: reduce_sum ((Uminus (Real a)) :: s)
       | (Integer a) :: b :: s -> b :: reduce_sum ((Integer a) :: s)
       | a :: s -> a :: reduce_sum s

  and simplify_times = fun
         (Real a) (Real b) -> (Real (a *. b))
       | (Integer a) (Integer b) -> (Integer (a * b))
       |  a (Real b) -> simplify_times (Real b) a 
       |  a (Integer b) -> simplify_times (Integer b) a
       |  (Uminus a) b -> simplify (Uminus (Times (a,b)))
       |  a (Uminus b) -> simplify (Uminus (Times (a,b)))
       | (Real a) (Times ((Real b), c)) -> 
	    simplify (Times ((Real (a *. b)), c))
       | (Integer a) (Times ((Integer b), c)) -> 
	    simplify (Times ((Integer (a * b)), c))
       | (Real a) b -> if (almost_equal a 0.0) then (Real 0.0)
                       else if (almost_equal a 1.0) then b
                       else if (almost_equal a (-1.0)) then simplify (Uminus b)
                       else Times ((Real a), b)
       | (Integer a) b -> if (a == 0) then (Integer 0)
                       else if (a == 1) then b
                       else if (a == (-1)) then simplify (Uminus b)
                       else Times ((Integer a), b)
       |  a b -> Times (a, b)
       
  and simplify_uminus = fun
         (Uminus x) -> x
       | (Plus a) -> 
            if (exists negative a)
	      then simplify (Plus (map negate a))
	    else Uminus (Plus a)
       | x -> Uminus x

  and simplify_real f = if (almost_equal f 0.0) then Real 0.0
             else if (f <. 0.0) then Uminus (Real (-. f))
	     else Real f

  and simplify_integer i = if (i < 0) then Uminus (Integer (- i))
	     else Integer i

  and negative = function
      Uminus x -> true
    | _ -> false

  and negate x = Uminus x

  and collect_minus x = 
      if (all negative x) then
         Uminus (Plus (map (fun (Uminus a) -> a) x))
      else
         Plus x

  and reorder = fun  (* push all Uminus to the end *)
      [] -> []
    | ((Uminus a) :: b) -> (reorder b) @ [(Uminus a)]
    | (a :: b) -> a :: (reorder b)

  and canonicalize = fun
      (Plus [a]) -> a
    | (Plus x) -> Plus (reorder x)
    | x -> x

  and filter_out f x = match x with 
	((Times (Real a, b)) as y :: c) -> 
		 let (q, r) = filter_out f c
		 in if (almost_equal a f) then (b::q, r)
		                else (q, y::r)
      | ((Uminus (Times (Real a, b))) as y :: c) -> 
		 let (q, r) = filter_out f c
		 in if (almost_equal a f) then ((Uminus b)::q, r)
		                else (q, y::r)
      |	_ -> ([], x)
		     

  and collect x = match x with 
      ((Times (Real a, b)) :: c) ->
		 let (witha, withouta) = filter_out a x
		 in (match witha with
		       [d] -> simplify (Times (Real a, d)) 
		              :: collect withouta
		     | _ -> simplify (Times (Real a, collect_minus witha))
		            :: collect withouta)
    | ((Uminus (Times (Real a, b))) :: c) ->
		 let (witha, withouta) = filter_out a x
		 in (match witha with
		       [d] -> simplify (Times (Real a, d)) :: collect withouta
		     | _ -> simplify (Times (Real a, collect_minus witha))
		            :: collect withouta)

    | (a :: c) -> a :: (collect c)
    | [] -> []

  and simplify_plus = fun
  	 [] -> Real 0.0  (* no terms, assume floating point *)
         | [a] -> a        (* one term *)
         | a -> canonicalize     (* everything else *)
              (collect_minus (collect a))

  in function  (* main simplify body *)
	 Var x -> Var x 
       | Subscript (a, b) -> Subscript (simplify a, simplify b)
       | Real f -> simplify_real f
       | Integer i -> simplify_integer i
       | Plus a ->  simplify_plus (reduce_sum (map simplify a))
       | Times (a, b) -> simplify_times (simplify a) (simplify b)
       | Comma (a, b) -> Comma (simplify a, simplify b)
       | Assign (a, b) -> Assign (simplify a, simplify b)
       | Call (a, b) -> Call (simplify a, simplify b)
       | Binop (op, a, b) -> Binop (op, simplify a, simplify b)
       | Uminus a -> simplify_uminus (simplify a)
;;

(*
 * cost model --- this simply counts the number of additions and
 * multiplications.  Returns a pair (additions, multiplications) 
 *)
let rec sum_costs = function
	 [] -> (0, 0)
       | ((a, b) :: s) -> let (c, d) = sum_costs s
                     in (a + c, b + d)
;;

let rec length = function
	 [] -> 0
       | _::s -> 1 + length s
;;

let rec cost_expr = function
         Var x -> (0, 0)
       | Subscript x -> (0, 0)
       | Real x -> (0, 0)
       | Integer x -> (0, 0)
       | Plus x -> let (a, m) = sum_costs (map cost_expr x)
	             in (a + (length x) - 1, m)
       | Times (a, b) -> let (a, m) = sum_costs [cost_expr a; cost_expr b]
	             in (a, m + 1)
       | Uminus a -> cost_expr a
       | Comma (a, b) -> sum_costs [cost_expr a; cost_expr b]
       | Assign (a, b) -> cost_expr b
       | Call (a, b) ->  sum_costs [cost_expr a; cost_expr b]
       | Binop a -> (0, 0)
;;

(*************************************
      program structure 
 *************************************)
type Declaration = Decl of string * Expr;; 

type Block = Block of (Declaration list) * (Stmt list)
and Stmt = Expr of Expr | Blk of Block
  | For of Expr * Expr * Expr * Block;;

type Function = Fcn of string * string * (Declaration list) * Block;;

let unparse_decl = function
    Decl (a, b) -> a ^ " " ^ unparse_expr b ^ ";\n"
;;

let foldr_string_concat = foldr (prefix ^) "";;

(* minor simplification: flatten nested blocks that do not declare variables *)
let rec simplify_block (Block (d, s)) = 
    let rec flatten_statement = function
	[] -> []
      | (Blk (Block ([], s))) :: q -> flatten_statement (s @ q)
      | (Blk (Block (a, s))) :: q ->
	     Blk (simplify_block (Block (a, s))) :: (flatten_statement q)
      | s :: q -> s :: (flatten_statement q)
   in
     Block (d, flatten_statement s)
;;

let id = "/* Generated by " ^ CVSId ^ " */\n\n";;

let rec unparse_block = function
  | Block (d, s) -> 
           "{\n"                                      ^ 
              foldr_string_concat (map unparse_decl d)    ^ 
              foldr_string_concat (map unparse_stmt s)    ^
           "}\n"
and unparse_stmt = function
         Expr x -> unparse_expr x ^ ";\n"
       | Blk b -> unparse_block b
       | For (a, b, c, d) ->
          "for (" ^
          unparse_expr a ^ "; " ^ unparse_expr b ^ "; " ^ unparse_expr c
          ^ ")" ^ unparse_block d
;;

let unparse_function = function
	 Fcn (typ, name, args, body) ->
	   let rec unparse_args = function
		[Decl (a, b)] -> a ^ " " ^ unparse_expr b 
	      | (Decl (a, b)) :: s -> a ^ " " ^ unparse_expr b  ^ ", "
				^  unparse_args s
	      | [] -> ""
	   in 
	       (typ ^ " " ^ name ^ "(" ^ unparse_args args ^ ")\n" ^
		 unparse_block body)
;;

let cost_decl x = (0, 0);;

let rec cost_block (Block (d, s)) = sum_costs (map cost_stmt s)
and cost_stmt = function
	 Expr x -> cost_expr x
       | Blk x -> cost_block x
       | For (a, b, c, x) -> cost_block x
;;

let cost_fcn (Fcn (a, b, c, x)) = cost_block x;;

(* canonicalization of constants *)
(*
 * The problem we are solving is the following: due to rounding
 * errors, numbers that should be the same are not.  For example,
 * sin(pi/4) and cos(pi/4).  In this case, the C compiler would hardwire
 * two constants instead of one, reducing speed.
 *
 * Therefore, we do one final pass over the expression tree
 * canonicalizing all the values (i.e., substituting almost-equal
 * numbers with equal)
 *)

let rec almost_mem k = function
     [] -> (false, k)
   | (k1 :: s) -> if (almost_equal k k1) then (true, k1)
                  else almost_mem k s
;;

(*
   canonicalize a constant and return the canonical value and the
   new set of constants
 *)
let canonicalize_konst konst konst_list =
    let (flag, k) = almost_mem konst konst_list
    in
	if flag then (k, konst_list)
	        else (konst, konst :: konst_list)
;;

let rec thread_list f kl l = match l with
      [] -> (kl, [])
    | (a :: s) -> 
         let (kl1, a1) = f kl a
	 in let (kl2, s1) = thread_list f kl1 s
	 in (kl2, a1 :: s1)
;;

let rec canonicalize_expr kl x = match x with
     Var a -> (kl, x)
    | Subscript a -> (kl, x)
    | Real a -> let (b, kl1) = canonicalize_konst a kl
	        in (kl1, Real b)
    | Integer a -> (kl, x)
    | Plus l -> let (kl1, l1) = thread_list canonicalize_expr kl l
	in (kl1, Plus l1)
    | Times (a, b) ->
          let (kl1, a1, b1) = thread kl a b
	  in (kl1, Times (a1, b1))
    | Uminus a -> let (kl1, a1) = canonicalize_expr kl a
	in (kl1, Uminus a1)
    | Comma (a, b) ->
          let (kl1, a1, b1) = thread kl a b
	  in (kl1, Comma (a1, b1))
    | Assign (a, b) ->
          let (kl1, a1, b1) = thread kl a b
	  in (kl1, Assign (a1, b1))
    | Call (a, b) ->
          let (kl1, a1, b1) = thread kl a b
	  in (kl1, Call (a1, b1))
    | Binop a -> (kl, x)

and thread kl a b = 
       let (kl1, a1) = canonicalize_expr kl a
       in let (kl2, b1) = canonicalize_expr kl1 b
       in (kl2, a1, b1)
;;

    
let rec canonicalize_block kl (Block (dl, sl)) =
    let (kl1, sl1) = thread_list canonicalize_stmt kl sl
    in (kl1, Block (dl, sl1))
and canonicalize_stmt kl = function
    Expr x -> let (kl1, x1) = canonicalize_expr kl x
	in (kl1, Expr x1)
  | Blk x -> let (kl1, x1) = canonicalize_block kl x
	in (kl1, Blk x1)
  | For (a, b, c, x) -> let (kl1, x1) = canonicalize_block kl x
        in (kl1, For (a, b, c, x1))
;;

let canonicalize_fcn (Fcn (a, b, c, d)) =
    let (k, d1) = canonicalize_block [] d
    in Fcn (a, b, c, d1)
;;

(****************************************
    complex numbers operations 
 ****************************************)
type Complex = Complex of Expr * Expr;;

let pi = 3.1415926535897932385;;

let COne = Complex (Real 1.0, Real 0.0);;

let CTimes (Complex (a,b)) (Complex (c,d)) = 
    Complex (Plus [Times (a,c); Uminus (Times (b,d))],
	     Plus [Times (a,d); Times (b,c)])
;;

let CUminus (Complex (a, b)) =
    Complex (Uminus a, Uminus b)
;;

(* complex exponential; returns exp(i * f) *)
let CExp f =
    Complex (Real (cos (2.0 *. pi *. f)), (Real (sin (2.0 *. pi *. f))))
;; 

(* complex assignment *)
let CAssign (Complex (ar, ai)) (Complex (br, bi)) =
    [Expr (Assign (ar, br)); Expr (Assign (ai, bi))]
;;

(* complex sum *)
let CPlus a =
    let rec unzip_complex = fun
    	[] -> ([], [])
  	| ((Complex (a,b)) :: s) ->
            let (r,i) = unzip_complex s
	  in
	    (a::r), (b::i)
    in
	let (c, d) = unzip_complex a
	in
	    Complex (Plus c, Plus d)
;;

(* complex simplification *)
let csimplify (Complex (a, b)) = Complex (simplify a, simplify b);;
let csimple (Complex (a, b)) = (simple a) && (simple b);;

(* extract real/imaginary *)
let c_real (Complex (a, b)) = Complex (a, Real 0.0);;
let c_imag (Complex (a, b)) = Complex (Real 0.0, b);;
let identity x = x;;

let c_conj (Complex (a, b)) = Complex (a, Uminus b);;

(* declaration of a complex variable *)
let CDecl (Complex (r, i)) =
    [Decl (real, r); Decl (real, i)]
;;

(* find the greatest factor of n smaller than sqrt(n) (1 if n is prime) *)
let factor n =
    let rec loop i f =
	if (i * i > n) then f
	else if ((n mod i) == 0) then loop (i + 1) i
	else loop (i + 1) f
    in
	loop 1 1
;;

(* fint the inverse of n modulo m *)
let invmod n m =
    let rec loop i =
	if ((i * n) mod m == 1) then i
	else loop (i + 1)
    in
	loop 1
;;

(* Yooklid's algorithm *)
let rec gcd n m =
    if (n > m)
      then gcd m n
    else
      let r = m mod n
      in
	  if (r == 0) then n
	  else gcd r n
;;
	    
(******************************************
          the fft generator 
 ******************************************)
let c_re = Var "c_re" ;;
let c_im = Var "c_im" ;;
let input = Var "in" ;;
let output = Var "out" ;;
let inout = Var "inout" ;;
let istride = Var "istride" ;;
let ostride = Var "ostride" ;;
let iostride = Var "stride" ;;
let W = Var "W";;

(* temporary variables *)
let gen_tmp depth n1 n2 =
    let tail = (string_of_int depth) ^ "_" ^ (string_of_int n1)
	^ "_" ^ (string_of_int n2)
    in let re = Var ("tre" ^ tail)
    and im = Var ("tim" ^ tail)
    in
	Complex (re, im)
;;

let declare_variables depth n1 n2 =
    forall_flat 0 n1 (fun j1 ->
	forall_flat 0 n2 (fun j2 ->
	    CDecl (gen_tmp depth j1 j2)))
;;


let rec fix_sign = fun
    (Complex (a, b)) (Complex (Uminus c, d)) ->
        fix_sign (Complex (Uminus a, b)) (Complex (c, d))
  | (Complex (a, b)) (Complex (c, Uminus d)) ->
        fix_sign (Complex (a, Uminus b)) (Complex (c, d))
  | a b -> (a, b)
;;

(* 
 * generate a temporary only if necessary (i.e., if expr is not simple)
 * Return a temporary declaration (maybe), code to assign the
 * temporary to the non-simple expression (maybe), and a simple expression.
 * If the non-simple expression has a - sign, assign the negative
 * to the temporary t and return -t.
 *)
let gen_tmp_maybe t expr =
    let sexpr = csimplify expr
    in
	if (csimple sexpr) then
	  ([], [], sexpr)
	else
	  let (t1, e1) = fix_sign t sexpr
	  in
	      (CDecl t, CAssign t (csimplify e1), t1)
;;

(* abstraction of sum_{i=0}^{n-1} *)
let Sigma a b f = csimplify (CPlus (forall_flat a b f));;

(* hack that transforms 
         Block ([Real variable, Imag variable...]
                [Real statement, Imag statement...])
   into 
         Block ([Real variable],
                [Real statement]),
         Block ([Imag variable],
                [Imag statement])

  (This turns out to be faster sometimes) 
 *)
let split_block (Block (d, s)) =
   let rec split_lists = fun
       [] -> [], []
       | (a :: b :: c) -> let (x, y) = split_lists c
                         in (a :: x), (b :: y)
   in let (d1, d2) = split_lists d
      and (s1, s2) = split_lists s
   in
      Block (d1, s1), Block (d2, s2)
;;

(*
 * Optimized FFTs for N = 3 and 5 (later, 7) from H. J. Nussbaumer, "Fast
 * Fourier Transform and Convolution Algorithms (2nd edition),"
 * (Springer-Verlag, New York, 1982). 
 * 
 * This code is disabled by default, since we have not
 * experimented with it enough yet...
 *)

let use_nussbaumer = false;;

let fftgen_3 input output sign depth =
    let newinput i = 
	let (d, a, e) = input i (gen_tmp depth i 0)
	  in e
    in let decls = forall_flat 0 3 (fun i ->
	let (d, a, e) = (input i (gen_tmp depth i 0))
          in d)
    and assignments = forall_flat 0 3 (fun i ->
	let (d, a, e) = (input i (gen_tmp depth i 0))
	  in a)
    and u = -1.0 *. sign *. 2.0 *. pi /. 3.0 
    in let (t1d, t1a, t1e) = (gen_tmp_maybe (gen_tmp (depth + 1) 1 0)
			      (CPlus [newinput 1; newinput 2]))
    in let m1 = CTimes (Complex (Real ((cos u) -. 1.0),Real 0.0)) t1e
    and (m0d, m0a, m0e) = (gen_tmp_maybe (gen_tmp (depth + 2) 0 0)
			   (csimplify (CPlus [newinput 0; t1e])))
    and (m2d, m2a, m2e) = (gen_tmp_maybe (gen_tmp (depth + 2) 2 0)
			   (csimplify
			    (CTimes (Complex (Real 0.0,Real (sin u)))
			     (CPlus [newinput 2; 
				     CUminus (newinput 1)]))))
    in let (s1d,s1a,s1e) = (gen_tmp_maybe (gen_tmp (depth + 3) 1 0)
			    (csimplify (CPlus [m0e;m1])))
    in let results = ((output 0 m0e) @ 
		      (output 1 (csimplify ((CPlus [s1e;m2e]))) @
    		      (output 2 (csimplify (CPlus [s1e; CUminus m2e])))))
    in Block (m0d @ s1d @ m2d,
	      [Blk (Block (decls,
			   assignments @ m2a @
			   [Blk (Block (t1d, t1a @ m0a @ s1a))]))] @
	      results)
;; 

let fftgen_5 input output sign depth =
    let newinput i = 
	let (d, a, e) = input i (gen_tmp depth i 0)
	  in e
    in let decls = forall_flat 0 5 (fun i ->
	let (d, a, e) = (input i (gen_tmp depth i 0))
          in d)
    and assignments = forall_flat 0 5 (fun i ->
	let (d, a, e) = (input i (gen_tmp depth i 0))
	  in a)
    and u = -1.0 *. sign *. 2.0 *. pi /. 5.0
    in let (t1d,t1a,t1e) = (gen_tmp_maybe (gen_tmp (depth + 1) 1 0)
			    (csimplify (CPlus [newinput 1;newinput 4])))
    and (t2d,t2a,t2e) = (gen_tmp_maybe (gen_tmp (depth + 1) 2 0)
			 (csimplify (CPlus [newinput 2;newinput 3])))
    and (t3d,t3a,t3e) = (gen_tmp_maybe (gen_tmp (depth + 1) 3 0)
			 (csimplify (CPlus [newinput 1;CUminus(newinput 4)])))
    and (t4d,t4a,t4e) = (gen_tmp_maybe (gen_tmp (depth + 1) 4 0)
			 (csimplify (CPlus [newinput 3;CUminus(newinput 2)])))
    in let (t5d,t5a,t5e) = (gen_tmp_maybe (gen_tmp (depth + 1) 5 0)
			    (csimplify (CPlus [t1e;t2e])))
    in let (m0d,m0a,m0e) = (gen_tmp_maybe (gen_tmp (depth + 2) 0 0)
			    (csimplify (CPlus [newinput 0;t5e])))
    and m1 = (CTimes (Complex (Real ((0.5 *. ((cos u) +. (cos (2.0 *. u))))
				     -. 1.0), Real 0.0))
	      t5e)
    and (m2d,m2a,m2e) = (gen_tmp_maybe (gen_tmp (depth + 2) 2 0)
			 (csimplify 
			  (CTimes (Complex (Real (0.5 *. ((cos u) -. 
							  (cos (2.0 *. u)))),
						     Real 0.0))
			   (CPlus [t1e;CUminus t2e]))))
    and (m3d,m3a,m3e) = (gen_tmp_maybe (gen_tmp (depth + 2) 3 0)
			 (csimplify
			  (CTimes (Complex (Real 0.0, (Uminus (Real (sin u)))))
			   (CPlus [t3e;t4e]))))
    and m4 = (CTimes (Complex (Real 0.0, 
			       (Real ((-1.0 *. (sin u)) -.
				      (sin (2.0 *. u))))))
	      t4e)
    and m5 = (CTimes (Complex (Real 0.0, 
			       (Real ((sin u) -.
				      (sin (2.0 *. u))))))
	      t3e)
    in let (s3d,s3a,s3e) = (gen_tmp_maybe (gen_tmp (depth + 3) 3 0)
			    (csimplify (CPlus [m3e;CUminus m4])))
    and (s5d,s5a,s5e) = (gen_tmp_maybe (gen_tmp (depth + 3) 5 0)
			 (csimplify (CPlus [m3e;m5])))
	
    and (s1d,s1a,s1e) = (gen_tmp_maybe (gen_tmp (depth + 3) 1 0)
			 (csimplify (CPlus [m0e;m1])))
    in let (s2d,s2a,s2e) = (gen_tmp_maybe (gen_tmp (depth + 3) 2 0)
			    (csimplify (CPlus [s1e;m2e])))	
    and (s4d,s4a,s4e) = (gen_tmp_maybe (gen_tmp (depth + 3) 4 0)
			 (csimplify (CPlus [s1e;CUminus m2e])))
    in let results = (
		      (output 1 (csimplify (CPlus [s2e;s3e]))) @
		      (output 2 (csimplify (CPlus [s4e;s5e]))) @
		      (output 3 (csimplify (CPlus [s4e;CUminus s5e]))) @
		      (output 4 (csimplify (CPlus [s2e;CUminus s3e]))))
    in (Block (s3d @ s5d @ s2d @ s4d,
	       [Blk (Block (decls,
			    assignments @
			    [Blk (Block (m3d @ t3d @ t4d,
					 t3a @ t4a @ m3a @
					 s3a @ s5a))] @
			    [Blk (Block (t5d @ m0d @ m2d,
					 [Blk (Block (t1d @ t2d,
						      t1a @ t2a @ t5a @ 
						      m2a))] @
				       m0a @ (output 0 m0e) @
				       [Blk (Block (s1d,
						    s1a @ s4a @ s2a))]))]))] @
	       results))
;;

(* generator for prime n: do n^2 algorithm for now (with some
   hackery for even/odd) *)
let fftgen_prime n input output sign depth =  
    let newinput i = 
	let (d, a, e) = input i (gen_tmp depth i 0)
	  in e
    in let decls = forall_flat 0 n (fun i ->
	let (d, a, e) = (input i (gen_tmp depth i 0))
          in d)
    and assignments = forall_flat 0 n (fun i ->
	let (d, a, e) = (input i (gen_tmp depth i 0))
	  in a)
    and sum i filter =
	Sigma 0 n (fun j ->
	    let fi = float_of_int i
	    and fj = float_of_int j
	    and fn = float_of_int n
	    in let coeff = filter (CExp(sign *. fi *. fj /. fn))
	    in [csimplify (CTimes coeff (newinput j))])
    in let computation_odd = 
	(output 0 (sum 0 identity)) @
	forall_flat 1 ((n + 1) / 2) (fun i ->
	    let (dr, ar, er) = gen_tmp_maybe (gen_tmp (depth + 1) 0 0)
                                             (sum i c_real)
	    and (di, ai, ei) = gen_tmp_maybe (gen_tmp (depth + 1) 1 0)
                                             (sum i c_imag)
	    in let (b1, b2) = 
                split_block
		    (Block (dr @ di,
			    ar @ ai @
			    (output i (CPlus [er; ei])) @
			    (output (n - i)
		             (csimplify (CPlus [er; CUminus ei])))))
            in [Blk b1; Blk b2])
    and computation_even =
	forall_flat 0 n (fun i ->
	    (output i (sum i identity)))
    in
	if (use_nussbaumer && n == 3) then
	  fftgen_3 input output sign depth
	else if (use_nussbaumer && n == 5) then
	  fftgen_5 input output sign depth
 	else if ((n mod 2) == 0) then
	  Block (decls, assignments @ computation_even)
	else
	  Block (decls, assignments @ computation_odd)
;;

(*
 * generator for the FFT blocks.
 *
 * n        : the size of the generated block
 * input    : int -> Expr, returns an expression containing the i-th input
 * output   : int -> Expr, returns an expression to which the i-th output
 *          : will be assigned
 *)
let rec fftgen n input output sign depth =
    let cooley_tukey n1 n2 =
	let decl = declare_variables depth n1 n2 
        and b1 = forall_flat 0 n2 (fun i2 ->
	    let newinput i1 =  input (i1 * n2 + i2)
	    and newoutput k1 = CAssign (gen_tmp depth k1 i2)
	    in  [Blk (fftgen n1 newinput newoutput sign (depth + 1))])
        and b2 = forall_flat 0 n1 (fun i1 ->
	    let fn = float_of_int n
	    and fi1 = float_of_int i1
	    in let newinput k2 ovar =  
		let fk2 = float_of_int k2
		in 
		    gen_tmp_maybe ovar 
   		         (CTimes (CExp (sign *. fi1 *. fk2 /. fn))
		                 (gen_tmp depth i1 k2)) 
		and newoutput k2 = output (i1 + k2 * n1)
	    in  [Blk (fftgen n2 newinput newoutput sign (depth + 1))])
      	in
	    Block (decl, b1 @ b2)

    and prime_factor n1 n2 =
	let t1 = invmod n1 n2
	and t2 = invmod n2 n1
	and decl = declare_variables depth n1 n2 
	in let b1 = forall_flat 0 n2 (fun i2 ->
	    let newinput i1 =  input ((i1 * n2 + i2 * n1) mod n)
	    and newoutput k1 = CAssign (gen_tmp depth k1 i2)
	    in  [Blk (fftgen n1 newinput newoutput sign (depth + 1))])
        and b2 = forall_flat 0 n1 (fun i1 ->
	    let newinput k2 ovar = ([], [], gen_tmp depth i1 k2)
	    and newoutput k2 = output ((i1 * t2 * n2 + k2 * t1 * n1) mod n)
	    in  [Blk (fftgen n2 newinput newoutput sign (depth + 1))])
      	in
	    Block (decl, b1 @ b2)

    in
   	let r = factor n
      	in
	    if (r == 1) (* n is prime *)
	    then 
	    	fftgen_prime n input output sign depth
	    else if ((gcd r (n / r)) == 1)
	    then 
	        prime_factor r (n/r)
	    else
	       	cooley_tukey r (n/r)
;;


(* generate the input/output variables *)
let gen_array_ref array stride k =
    let imult a b = Times ((Integer a), b)
    in
        csimplify (Complex (re, im))
    where re = Call(c_re, Subscript(array, imult k stride))
    and   im = Call(c_im, Subscript(array, imult k stride))
;;

let make_input_function array stride k ovar =
  (CDecl ovar,
  CAssign ovar (gen_array_ref array stride k),
  ovar)
;;

let make_output_function array stride k =
    CAssign (gen_array_ref array stride k)
;;

let gen_input = make_input_function input istride;;
let gen_output = make_output_function output ostride;;
let gen_inout_in = make_input_function inout iostride;;
let gen_inout_out = make_output_function inout iostride;;

type Direction = FORWARD | BACKWARD;;

let gen_inout_twiddle dir i ovar =
  let conj_maybe = match dir with
      FORWARD -> identity
    | BACKWARD -> c_conj
  in
  if (i == 0) then
       gen_inout_in i ovar 
  else
       let t = Complex (Var "tr", Var "ti")
       and tw = Complex (Var "twr", Var "twi")
       in let x = csimplify (CTimes t (conj_maybe tw))
       in let (ovar1, x1) = fix_sign ovar x
    in
    	(CDecl ovar,
         [Blk (Block (CDecl t @ CDecl tw, 
                CAssign t (gen_array_ref inout iostride i) @
                CAssign tw (gen_array_ref W (Integer 1) (i - 1)) @
                CAssign ovar x1))],
         ovar1)
;;

let gen_block n input output sign =
    simplify_block (fftgen n input output sign 0)
;;

let print_cost f =
    let (a, m) = cost_fcn f
    in
	"/* This function contains " ^
	(string_of_int a) ^ " FP additions and "  ^
	(string_of_int m) ^ " FP multiplications */\n\n"
;;

let no_twiddle_gen n dir =
    let ns = string_of_int n
    and (name, sign) = match dir with
	  FORWARD -> "fftw_no_twiddle_", (-. 1.0)
	| BACKWARD -> "fftwi_no_twiddle_", 1.0
    in
     let tree = 
	Fcn ("void", name ^ ns,
	     [Decl ("const FFTW_COMPLEX *", input);
	      Decl ("FFTW_COMPLEX *", output);
	      Decl ("int", istride);
	      Decl ("int", ostride)],
	     (gen_block n gen_input gen_output sign))
     in id ^ (print_cost tree) ^ (unparse_function (canonicalize_fcn tree))
;;

let twiddle_gen n dir =
    let ns = string_of_int n
    and m = Var "m";
    and dist = Var "dist";
    and A = Var "A"
    and i = Var "i"
    and (name, sign) = match dir with
	  FORWARD -> "fftw_twiddle_", (-. 1.0)
	| BACKWARD -> "fftwi_twiddle_", 1.0
    in
      let body = Block (
	    [Decl ("int", i);
	     Decl ("COMPLEX *", inout)],
	    [Expr (Assign (inout, A));
             For (Assign (i, Integer 0),
		  Binop (" < ", i, m),
		  Comma (Assign (i, Plus [i; Integer 1]),
			(Comma (Assign (inout, Plus [inout; dist]),
                               (Assign (W, Plus [W; Integer (n - 1)]))))),
                  (gen_block n (gen_inout_twiddle dir) gen_inout_out sign))])

     in let tree = 
	Fcn ("void", name ^ ns,
	     [Decl ("FFTW_COMPLEX *", A);
	      Decl ("const FFTW_COMPLEX *", W);
	      Decl ("int", iostride);
	      Decl ("int", m);
	      Decl ("int", dist)],
              body)

     in id ^ (print_cost tree) ^ (unparse_function (canonicalize_fcn tree))
;;

(**********************************
      main program
 **********************************)
let usage = "Usage: genfft [-notwiddle | -notwiddleinv | -twiddle | -twiddleinv] <number>";;

if sys__interactive then () else
if vect_length sys__command_line <> 3 then begin
  print usage;
  print_newline()
end else begin
  try
      if sys__command_line.(1) = "-notwiddle" then
	print (no_twiddle_gen (int_of_string sys__command_line.(2)) FORWARD)
      else
      if sys__command_line.(1) = "-notwiddleinv" then
	print (no_twiddle_gen (int_of_string sys__command_line.(2)) BACKWARD)
      else
      if sys__command_line.(1) = "-twiddle" then
	print (twiddle_gen (int_of_string sys__command_line.(2)) FORWARD)
      else
      if sys__command_line.(1) = "-twiddleinv" then
	print (twiddle_gen (int_of_string sys__command_line.(2)) BACKWARD)
      else
	print_string usage;
    print_newline()
  with Failure "int_of_string" ->
    print_string "Bad integer constant";
    print_newline()
end
;;
