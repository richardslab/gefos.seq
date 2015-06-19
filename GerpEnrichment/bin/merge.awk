BEGIN {
    while(( getline line < fn) > 0 ) {
	split(line,items,"\t")
	pvalues[items[3]] = items[12]
    }
    print length(pvalues) > "/dev/stderr"
}

{
    if ($6 in pvalues){
	print $0, pvalues[$6]
    }
}
