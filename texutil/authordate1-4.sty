%  This LaTeX style-option file, authordate1-4.sty, is intended for documents 
%  that use the author-date citation system.  When used in conjunction with
%  authordate1.bst, ... , authordate4.bst, it allows citations of the
%  form \shortcite{bloggs-60} as well as the usual \cite{bloggs-60}.
%  When
%       ... \cite{bloggs-60} ... \shortcite{bloggs-60} ...
%  appears in the input file,
%       ... (Bloggs, 1960) ... (1960) ...
%  appears in the final document.
%
%  Note that:
%  - a command \bibtitle is defined, to generate the bibliography title
%  - a command \bibheadtitle is defined, to generate the text to be used
%    by \pagestyle{headings}
%  - \bibtitle and \bibheadtitle are set to generate "References" and
%    "REFERENCES" respectively.
%  - an \addcontentsline{toc}{...}{\bibtitle} is issued.
%  Thus:
%  - the default bibliography title is "References", but can be changed by
%    having \renewcommand{\bibtitle}{...} before \begin{thebibliography}
%    or \bibliography
%  - the default page header text is "REFERENCES", but can be changed by
%    \renewcommand{\bibheadtitle}{...}                  
%  - \tableofcontents will give a table of contents that mentions the
%    bibliography.
%
%                                                     David Rhead
%                                                     Cripps Computing Centre
%                                                     February 1990

%  Code taken from apalike.sty.
%
\def\@cite#1#2{(#1\if@tempswa , #2\fi)}
\def\@biblabel#1{}
 
% Set length of hanging indentation for bibliography entries
%
\newlength{\bibhang}
\setlength{\bibhang}{2em}

%  Define default title and page-head text.
%
\def\bibtitle{References}
\def\bibheadtitle{REFERENCES}
 
\@ifundefined{chapter}{\def\thebibliography#1{\section*{\bibtitle\@mkboth
  {\bibheadtitle}{\bibheadtitle}}
  \addcontentsline{toc}{section}{\bibtitle}\list            %  For styles in
  {\relax}{\setlength{\labelsep}{0em}                       %  which \chapter
        \setlength{\itemindent}{-\bibhang}                  %  is undefined.
        \setlength{\leftmargin}{\bibhang}}
    \def\newblock{\hskip .11em plus .33em minus .07em}
    \sloppy\clubpenalty4000\widowpenalty4000
    \sfcode`\.=1000\relax}}%
{\def\thebibliography#1{\chapter*{\bibtitle\@mkboth           %  For styles in
  {\bibheadtitle}{\bibheadtitle}}                             %  which \chapter
  \addcontentsline{toc}{chapter}{\bibtitle}\list              %  is defined.
  {\relax}{\setlength{\labelsep}{0em}
        \setlength{\itemindent}{-\bibhang}
        \setlength{\leftmargin}{\bibhang}}
    \def\newblock{\hskip .11em plus .33em minus .07em}
    \sloppy\clubpenalty4000\widowpenalty4000
    \sfcode`\.=1000\relax}}
%
%  Code taken from aaai.sty.
%
% don't box citations, add space between multiple citations, separate with ;
\def\@citex[#1]#2{\if@filesw\immediate\write\@auxout{\string\citation{#2}}\fi
  \def\@citea{}\@cite{\@for\@citeb:=#2\do
    {\@citea\def\@citea{; }\@ifundefined
       {b@\@citeb}{{\bf ?}\@warning
       {Citation `\@citeb' on page \thepage \space undefined}}%
{\csname b@\@citeb\endcsname}}}{#1}}
% Allow short (name-less) citations, when used in
% conjunction with a bibliography style that creates labels like
%       \citename{<names>, }<year>
\let\@internalcite\cite
\def\cite{\def\citename##1{##1}\@internalcite}
\def\shortcite{\def\citename##1{}\@internalcite}
