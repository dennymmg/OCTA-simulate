function [Size]=getParticleSize(depth1,df,nbp,Size_p)
    if (isPartiallyLeft(depth1,df,nbp,Size_p))
        Size=nbp+Size_p-depth1                // contribution of size of particle to bin

    elseif(isPartiallyRight(depth1,df,nbp,Size_p))
        Size=(depth1+df)-nbp                   //contribution of size of particle to bin

    elseif(isCompletelyInside(depth1,df,nbp,Size_p))
        Size=Size_p                            //contribution of size of particle to bin

    elseif(isCompletelyRightOutside(depth1,df,nbp,Size_p))
        Size=0                              //since the particle is completely outside the bin 
                        //thus the contribution of size  to bin is zero 
    end
endfunction

function [z]=isPartiallyLeft(depth1,df,nbp,Size)
    if (depth1>nbp && nbp+Size>depth1)
        z = %T
    else
        z = %F
    end
endfunction

function [z]=isPartiallyRight(depth1,df,nbp,Size)
    if (depth1+df>nbp && nbp+Size>depth1+df)
        z = %T
    else
        z = %F
    end
endfunction

function [z]=isCompletelyInside(depth1,df,nbp,Size)
    if (depth1<=nbp && nbp+Size<=depth1+df)
        z = %T
    else
        z = %F
    end
endfunction

function [z]=isCompletelyRightOutside(depth1,df,nbp,Size)
    if (depth1+df<nbp )
        z = %T
    else
        z = %F
    end
endfunction

