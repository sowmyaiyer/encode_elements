{
        split($6, intervalIds, ",")
        split($5, scores, ",")
        num=length(intervalIds)
        if (num == 1)
        {
		split(intervalIds[1],interval,"_")
                print($1"\t"$2"\t"$3"\t"$7"\t"$5"\t"$4)
        }
        else {
                argmax=-1
                max=-1
                for (i = 1 ; i <= num; i ++)
                {
                        if (scores[i] > max)
                        {
                                max=scores[i]
                                argmax=i
                        }
                }
                split(intervalIds[argmax], interval,"_")
                print(interval[1]"\t"interval[2]"\t"interval[3]"\t"$7"\t"max"\t"$4)
        }
}
