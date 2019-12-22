\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{mrm2plain}[2018/07/21, v1.0]

\LoadClass[a4paper]{article}
%\LoadClass{scrartcl}
\RequirePackage[numbers]{natbib}
\RequirePackage{graphicx}
\RequirePackage{amsmath}
\RequirePackage{mathtools}
\RequirePackage{amssymb}

% font stuff
\RequirePackage{microtype}
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}
% \RequirePackage{bera}

\newcommand{\papertype}[1]{\def\@papertype{#1}}
\newcommand{\corraddress}[1]{\def\@corraddress{#1}}
\newcommand{\corremail}[1]{\def\@corremail{#1}}
\newcommand{\keywords}[1]{}
\newcommand{\runningauthor}[1]{}

\RequirePackage{authblk}
\RequirePackage{xpatch}
\patchcmd{\AB@affilsepx}{\\}{\\[3\p@]}{}{}
\patchcmd{\@author}{\AB@authlist\\[\affilsep]\AB@affillist}{\AB@authlist}{}{}

\def\@maketitle{%
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
  \let \footnote \thanks
    {\LARGE \@title \par}%
    \vskip 1.5em%
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
    %\vskip 1em%
    %{\large \@date}%
  \end{center}%
  \par
  \AB@affillist\par%
  \vskip 1.5em}

