BEGIN {
    while(( getline line < fn) > 0 ) {
	split(line,items," ")
	snps[items[2]] = 1
    }
    print length(snps) > "/dev/stderr"
}

{
    if (NR == 1) { print $0 }
    if ($6 in snps){
	print $0
    }
}