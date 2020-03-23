function index=x(i,j,f)
    global l_max
    global f_total
    index=(i-1)*l_max*f_total+(j-1)*f_total+f;
end