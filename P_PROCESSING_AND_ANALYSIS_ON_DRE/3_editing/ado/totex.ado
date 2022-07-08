*! version 0.9  9Mar2006 rosa gini
*! version 1.0  20jan2014 rosa gini

program define totex,rclass
	version 8.2
	args strnotex
	*local strnotex: subinstr local strnotex "\" "$\backslash$",all
	local strnotex: subinstr local strnotex "%" "\%",all
	local strnotex: subinstr local strnotex "$" "\$",all
	local strnotex: subinstr local strnotex "_" "\_",all
	local strnotex: subinstr local strnotex "&" "\&",all
	local strnotex: subinstr local strnotex "#" "\#",all
	return local stringtex `"`strnotex'"'
end

