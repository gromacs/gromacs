%% This is file `rugletter.cls', adapted from
%% `letter.cls' 1995/07/06 v1.2s by Ronald Tol and Rein Smedinga
%% (C) 1995 vakgroep Informatica RUG.
%%
%% This class automatically included RuG2e.sty and Babel which are
%% needed for some at least some of the letters of the department.
%%
%% 1996/06/24 updates interaction with RuG2e.sty (for eps-option)
%% 1996/04/03 replaces epsf-package by graphics-package to get a standard
%% 1996/03/20 created RuG2eletter.cls
%% 1996/05/08 renamed the class to rugletter.cls
%%
\NeedsTeXFormat{LaTeX2e}[1995/06/01]
\ProvidesClass{rugletter}
              [1996/03/20 v1.1
               RuG LaTeX document class]
\newcommand\@ptsize{}
% only A4 possible in RuG2eletter
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\if@compatibility
  \DeclareOption{twoside}{\@latexerr{No `twoside' layout for letters}%
                                   \@eha}
\else
  \DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\fi
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{eps}{\PassOptionsToPackage{\CurrentOption}{RuG2e}}
\DeclareOption*{\PassOptionsToPackage{\CurrentOption}{babel}}
\ExecuteOptions{a4paper,11pt,oneside,onecolumn,final}
\ProcessOptions
\RequirePackage{babel}
\RequirePackage{graphics}
\input{size1\@ptsize.clo}
\RequirePackage{RuG2e}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0.7em}
\setlength\parindent{0\p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setlength\headheight{12\p@}
\setlength\headsep   {45\p@}
\setlength\footskip{25\p@}
\if@compatibility
  \setlength\textwidth{365\p@}
  \setlength\textheight{505\p@}
\fi
\if@compatibility
  \setlength\oddsidemargin{53pt}
  \setlength\evensidemargin{53pt}
  \setlength\marginparwidth{90pt}
\else
  \setlength\@tempdima{\paperwidth}
  \addtolength\@tempdima{-2in}
  \addtolength\@tempdima{-\textwidth}
  \setlength\oddsidemargin   {.5\@tempdima}
  \setlength\evensidemargin  {\oddsidemargin}
  \setlength\marginparwidth  {90\p@}
\fi
\setlength\marginparsep {11\p@}
\setlength\marginparpush{5\p@}
\setlength\topmargin{0pt}
\setlength\footnotesep{12\p@}
\setlength{\skip\footins}{10\p@ \@plus 2\p@ \@minus 4\p@}


%%%%%%%%%%%%%%%%%%%% headings
\def\ps@headings{\iflanguage{dutch}%
{\def\@oddhead{\iftelefbrief \hfil\@telefaxinfo\hfil
                                \else \small\hifnt%
\headtoname\ \ignorespaces\toname \hfil \@date \hfil \pagename\ \thepage\fi}\def\@oddfoot{}}
{\def\@oddhead{\iftelefbrief \hfil\@telefaxinfo\hfil 
                                \else \small\hifnt%
\headtoname\ \ignorespaces\toname \hfil \@date \hfil \pagename\ \thepage\fi}\def\@oddfoot{}}
}
\def\ps@empty{%
      \let\@oddfoot\@empty\let\@oddhead\@empty
      \let\@evenfoot\@empty\let\@evenhead\@empty}
\def\ps@firstpage{%
      \def\@oddhead{\letterhead\hfill}\def\@oddfoot{}
      \def\@evenhead{}\def\@evenfoot{}}
\def\ps@plain{%
      \let\@oddhead\@empty
      \def\@oddfoot{\normalfont\hfil\thepage\hfil}%
      \def\@evenfoot{\normalfont\hfil\thepage\hfil}}
\newcommand*{\name}[1]{\def\fromname{#1}}
\newcommand*{\signature}[1]{\def\fromsig{#1}}
\newcommand*{\address}[1]{\def\fromaddress{#1}}
\newcommand*{\location}[1]{\def\fromlocation{#1}}
\newcommand*{\telephone}[1]{\def\telephonenum{#1}}
\newcommand*{\email}[1]{\def\emailname{#1}}
\newcommand*{\ccname}{cc}
\newcommand*{\enclname}{encl}
\newcommand*{\pagename}{Page}

\name{}
\signature{}
\address{}
\location{}
\telephone{}
\email{}

%%%%%%%%%%%% RuG-letter dependent macros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% definition of the name of the place (normally only in Dutch letters)
\def\Placename{\iflanguage{dutch}%
{Groningen, }
{\relax}
}
\def\plaats#1{\def\Placename{#1, }}

% definition of the date-line
\def\@datumregel{\Placename\@date.}

\def\Betreftname{\iflanguage{dutch}%
{Betreft}
{Subject}
}

\def\@Ref{\iflanguage{dutch}%
{Ref.}
{Ref.}
}

%speciaal voor Nederlandse brieven:
\def\plaats#1{\def\@Place{#1, }}
\def\bijlage#1{\par\noindent
\parbox[t]{\textwidth}{\@hangfrom{\rm Bijlage: }\ignorespaces #1\strut}\par}
\def\bijlagen#1{\par\noindent
\parbox[t]{\textwidth}{\@hangfrom{\rm Bijlagen: }\ignorespaces #1\strut}\par}
\def\betreft#1{\def\@betreft{#1}}
\def\@betreft{}
\def\reference#1{\def\@reference{#1}}
\def\@reference{}
\def\owninfo#1{\def\@owninfo{#1}}
\let\subject\betreft

\def\makelabels{%
  \AtBeginDocument{%
     \let\@startlabels\startlabels
     \let\@mlabel\mlabel
     \if@filesw
       \immediate\write\@mainaux{\string\@startlabels}\fi}%
  \AtEndDocument{%
     \if@filesw\immediate\write\@mainaux{\string\clearpage}\fi}}
\@onlypreamble\makelabels

\newcount\letno \letno=0
\newenvironment{letter}[1]
  {\newpage
   \global\advance\letno by 1
    \@ifundefined{let-\the\letno}%
    {\global\def\totalpages{1}}%
     {\global\def\totalpages{\expandafter\csname let-\the\letno\endcsname}} 
    \if@twoside \ifodd\c@page
                \else\thispagestyle{empty} \hbox{}\newpage\fi
    \fi
    \c@page\@ne
    \interlinepenalty=200 % smaller than the TeXbook value
    \@processto{\leavevmode\ignorespaces #1}}
  {\stopletter
    \if@filesw
      \begingroup
        \let\\=\relax
        \let\protect\@unexpandable@protect
        \immediate\write\@auxout{%
          \string\@mlabel{\returnaddress}{\toname\\\toaddress}}%
        \immediate\write\@auxout{%
           \string\global\string\expandafter\string\def\string\csname\space
           let-\the\letno\string\endcsname{\thepage}}
           \if\thepage\totalpages \fi 
      \endgroup
    \fi
    \@@par\pagebreak\@@par}
\long\def\@processto#1{\@xproc #1\\@@@\ifx\toaddress\@empty
    \else \@yproc #1@@@\fi}
\long\def\@xproc #1\\#2@@@{\def\toname{#1}\def\toaddress{#2}}
\long\def\@yproc #1\\#2@@@{\def\toaddress{#2}}
\newcommand*{\stopbreaks}{\interlinepenalty \@M
   \def\par{\@@par\nobreak}\let\\=\@nobreakcr
   \let\vspace\@nobreakvspace}
\def\@nobreakvspace{\@ifstar{\@nobreakvspacex}{\@nobreakvspacex}}

\def\@nobreakvspacex#1{\ifvmode\nobreak\vskip #1\relax\else
               \@bsphack\vadjust{\nobreak\vskip #1}\@esphack\fi}

\def\@nobreakcr{\vadjust{\penalty\@M}\@ifstar{\@xnewline}{\@xnewline}}
\newcommand*{\startbreaks}{\let\\=\@normalcr
   \interlinepenalty 200\def\par{\@@par\penalty 200\relax}}

%%%%%%%%%%%%%%%%% telefax %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newif\iftelefbrief
\global\telefbrieffalse
%
\def\telef{\global\telefbrieftrue\pagestyle{headings}
\large
}
\def\endtelef{\global\telefbrieffalse}

\def\@telefaxinfo{TELEFAX -- \pagename\ \thepage\ \ofname\ \totalpages --}

\def\ofname{\iflanguage{dutch}%
{van}
{of}
}

\def\fax{\@ifnextchar [ {\faxi}{\faxii}}
\def\faxi[#1]#2{\telef\letter{TELEFAX #1 \\\strut\\ #2}}
\def\faxii{\telef\letter}
\def\endfax{\endletter\endtelef}


\newdimen\longindentation
\longindentation=.5\textwidth
\newdimen\indentedwidth
\indentedwidth=\textwidth
\newcommand*{\opening}[1]{\ifx\@empty\fromaddress
    \thispagestyle{firstpage}%
    \vspace*{-7mm}% changed from 1mm  RS 24-8-93
    \llap{%%%address box moved to left because of bigger window in envelope
%%%   RS 15-12-94
    \hbox to 2cm{\vbox to 45mm{.\strut\vfil\noindent
    \begin{tabular}{l}
    \toname \\
    \toaddress
    \end{tabular}\vfil.\hfill
    }%end vbox
    \hss}%end hbox
    }%% end llap
    \vspace*{10mm}% changed from 3mm  RS 24-8-93
    \strut\hfill \@datumregel%
    \par\vspace{2\parskip}%\noindent%
  \else  % home address
    \thispagestyle{empty}%
    \vspace*{-10\parskip}
    {\raggedright\ignorespaces\fromaddress }\\*[8\parskip]
    {\raggedright\ignorespaces
    \toname \\ \toaddress 
    }\par\vspace{2\parskip}
    {\strut\hfill \@datumregel}%
    \par\vspace{5\parskip}
  \fi
  \ifx\@betreft\@empty\else
  \parbox[t]{\textwidth}{\@hangfrom{\rm \Betreftname: }\ignorespaces \@betreft}%  % changed in reference (RMT)
  \ifx\@reference\@empty\par\vspace{2\parskip}\else\\\fi
  \fi
  \ifx\@reference\@empty\else
  \parbox[t]{\textwidth}{\@hangfrom{\rm \@Ref: }\ignorespaces \@reference}%
  \par\vspace{2\parskip}\fi
  \par#1\par\nobreak}

\newcommand{\closing}[1]{\par\nobreak\vspace{\parskip}%
  \stopbreaks
  \vspace{2\parskip}
  \noindent
  \parbox{\longindentation}{\raggedright
  \ignorespaces #1\strut\\%%%%\strut added 8-8-90
  \vspace{4\baselineskip}%%%%4\baselineskip instead of 8\parskip
  \ifx\@empty\fromsig
  \fromname
  \else \fromsig \fi\strut}
  \par\vfill}

\medskipamount=\parskip
\newcommand*{\cc}[1]{%
  \par\noindent
  \parbox[t]{\textwidth}{%
    \@hangfrom{\normalfont\ccname: }%
    \ignorespaces #1\strut}\par}
\newcommand*{\encl}[1]{%
  \par\noindent
  \parbox[t]{\textwidth}{%
    \@hangfrom{\normalfont\enclname: }%
    \ignorespaces #1\strut}\par}
\newcommand*{\ps}{\par\startbreaks}
\newcommand*{\stopletter}{}
\newcommand*{\returnaddress}{}
\newcount\labelcount
\newcommand*{\startlabels}{\labelcount\z@
  \pagestyle{empty}%
  \let\@texttop\relax
  \topmargin -50\p@
  \headsep \z@
  \oddsidemargin -35\p@
  \evensidemargin -35\p@
  \textheight 10in
  \@colht\textheight  \@colroom\textheight \vsize\textheight
  \textwidth 550\p@
  \columnsep 26\p@
  \ifcase \@ptsize\relax
    \normalsize
  \or
    \small
  \or
    \footnotesize
  \fi
  \baselineskip \z@
  \lineskip \z@
  \boxmaxdepth \z@
  \parindent \z@
  \twocolumn\relax}
\let\@startlabels=\relax
\newcommand*{\mlabel}[2]{%
  \parbox[b][2in][c]{262\p@}{\strut\ignorespaces #2}%
  }
\let\@mlabel=\@gobbletwo
\setlength\leftmargini  {2.5em}
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\setlength\leftmarginv  {1em}
\setlength\leftmarginvi {1em}
\setlength\leftmargin    {\leftmargini}
\setlength  \labelsep  {5\p@}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\setlength\partopsep{0\p@}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\def\@listI{\setlength\leftmargin{\leftmargini}
            \setlength\parsep {0\p@}%
            \setlength\topsep {.4em}%
            \setlength\itemsep{.4em}}
\let\@listi\@listI
\@listi
\def\@listii {\setlength  \leftmargin{\leftmarginii}%
              \setlength  \labelwidth{\leftmarginii}%
              \addtolength\labelwidth{-\labelsep}}
\def\@listiii{\setlength  \leftmargin{\leftmarginiii}%
              \setlength  \labelwidth{\leftmarginiii}%
              \addtolength\labelwidth{-\labelsep}%
              \setlength  \topsep    {.2em}%
              \setlength  \itemsep   {\topsep}}
\def\@listiv {\setlength  \leftmargin{\leftmarginiv}%
              \setlength  \labelwidth{\leftmarginiv}%
              \addtolength\labelwidth{-\labelsep}}
\def\@listv  {\setlength  \leftmargin{\leftmarginv}%
              \setlength  \labelwidth{\leftmarginv}%
              \addtolength\labelwidth{-\labelsep}}
\def\@listvi {\setlength  \leftmargin{\leftmarginvi}%
              \setlength  \labelwidth{\leftmarginvi}%
              \addtolength\labelwidth{-\labelsep}}
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{$\m@th\bullet$}
\newcommand\labelitemii{\normalfont\bfseries --}
\newcommand\labelitemiii{$\m@th\ast$}
\newcommand\labelitemiv{$\m@th\cdot$}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*{\descriptionlabel}[1]{\hspace\labelsep
                                \normalfont\bfseries #1}

\newenvironment{verse}
               {\let\\=\@centercr
                \list{}{\setlength\itemsep{\z@}%
                        \setlength\itemindent{-15\p@}%
                        \setlength\listparindent{\itemindent}%
                        \setlength\rightmargin{\leftmargin}%
                        \addtolength\leftmargin{15\p@}}%
                \item[]}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\setlength\listparindent{1.5em}%
                        \setlength\itemindent{\listparindent}%
                        \setlength\rightmargin{\leftmargin}}%
                \item[]}
               {\endlist}
\newenvironment{quote}
               {\list{}{\setlength\rightmargin{\leftmargin}}%
                \item[]}
               {\endlist}

% bibliography-environment also possible in letter environment
\def\thebibliography#1{\par{\bf \bibname}
\list{[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}\leftmargin\labelwidth
 \advance\leftmargin\labelsep
 \usecounter{enumi}}
 \def\newblock{\hskip .11em plus .33em minus -.07em}
 \sloppy
 \sfcode`\.=1000\relax}
\let\endthebibliography=\endlist
\def\bibname{Bibliography}

\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand\theequation{\@arabic\c@equation}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\relax}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\relax}
\DeclareRobustCommand*{\cal}{\@fontswitch{\relax}{\mathcal}}
\DeclareRobustCommand*{\mit}{\@fontswitch{\relax}{\mathnormal}}
\renewcommand\footnoterule{%
  \kern-\p@
  \hrule \@width .4\columnwidth
  \kern .6\p@}
\long\def\@makefntext#1{%
    \noindent
    \hangindent 5\p@
    \hb@xt@5\p@{\hss\@makefnmark}#1}


\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{plain}
\pagenumbering{arabic}
\raggedbottom
\let\@texttop=\relax
\onecolumn
\endinput
%% 
%% End of file `RuG2eletter.cls'.


