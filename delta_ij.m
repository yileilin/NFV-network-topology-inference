function index=delta_ij(f1,f2)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+(f1-1)*f_total+f2;
end