function p_k = regular_wf( g_k,lamda )
%Regular WF
if lamda<g_k
    p_k=1/lamda-1/g_k;
else
    p_k=0;
end
end

