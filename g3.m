function index=g3(f1,f2,A)
    global n_path
    global l_max
    global f_total
    global n_dummy
    global n_category
    index=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        f_total*f_total*l_max*l_max*n_path+f_total*f_total*n_path+ ...,
        (f1-1)*f_total*n_category+(f2-1)*n_category+A;
end