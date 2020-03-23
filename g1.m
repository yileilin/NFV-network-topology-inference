function index=g1(f1,f2,j1,j2,i)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        (f1-1)*f_total*l_max*l_max*n_path+(f2-1)*l_max*l_max*n_path+ ...,
        (j1-1)*l_max*n_path+(j2-1)*n_path+i;
end