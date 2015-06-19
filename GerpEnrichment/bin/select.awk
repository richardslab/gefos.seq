BEGIN {
    prev_tag=""
    min_p=1
    min_p_tag=""
}

{
    tag=$3
    if (NR == 1){
	min_p_tag=tag
	min_p=$8
    }
    if ((prev_tag != "") && (tag != prev_tag)){
	print prev_tag, min_p_tag, min_p
	min_p_tag=$6
	min_p=$8
    }else{
	if (tag == prev_tag){
	    if ($8 < min_p){
		min_p=$8
		min_p_tag=$6
	    }
	}else{
	    printf "ERROR: tag:'%s' prev_tag:'%s'\n", tag, prev_tag > "/dev/stderr"
	}
    }
    prev_tag=tag
}