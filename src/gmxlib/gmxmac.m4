ifdef(`USEF77',`define(`SCALARG',`*$1')',`define(`SCALARG',`$1')')
ifdef(`USEF77',`define(`SCAL',`&($1)')',`define(`SCAL',`$1')')
ifdef(`USEF77',`define(`FUNC',FUNCTION(`$1'))',`define(`FUNC',`$1')')

