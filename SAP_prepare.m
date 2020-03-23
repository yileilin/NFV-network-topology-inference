%define some variables
global n_path
global l_max
global n_dest
global f_total
global n_dummy
global n_source

N_VNF=5;
n_path=6;
l_max=12;
n_dummy=5;
n_dest=3;
n_source=2;
f_total=n_source+n_dest+N_VNF+n_dummy; %f_total=1+n_dest+N_VNF+n_dummy

Aneq=sparse(1,100);
Bneq=sparse(1,1);
i_neq_left=1;
%(8)
for i=1:1:n_path
    for j=1:1:l_max
        Aneq(i_neq_left,x(i,j,1):x(i,j,f_total))=1;
        i_neq_left=i_neq_left+1;
    end
end
Bneq(1:(n_path*l_max),1)=1;
i_neq_right=n_path*l_max+1;
i_neq_left=i_neq_right;
        
%(13)
for f=n_source+1+n_dest+N_VNF:1:f_total
    for i=1:1:n_path
        for j=1:1:l_max
            Aneq(i_neq_left,x(i,j,f))=1;
            Aneq(i_neq_left,delta_k(f-n_source-n_dest-N_VNF))=-1;
            i_neq_left=i_neq_left+1;
        end
    end
end

Bneq(i_neq_right:(i_neq_left-1),1)=0;
i_neq_right=i_neq_right+1;

%(15)
for i=1:1:n_path
    for f1=1:1:f_total
        for f2=1:1:f_total
            for j1=1:1:l_max
                for j2=1:1:l_max
                    Aneq(i_neq_left,g1(f1,f2,j1,j2,i))=1;
                    Aneq(i_neq_left,x(i,j1,f1))=-1;
                    i_neq_left=i_neq_left+1;
                    Aneq(i_neq_left,g1(f1,f2,j1,j2,i))=1;
                    Aneq(i_neq_left,x(i,j2,f2))=-1;
                    i_neq_left=i_neq_left+1;
                    Aneq(i_neq_left,g1(f1,f2,j1,j2,i))=-1;
                    Aneq(i_neq_left,x(i,j1,f1))=1;
                    Aneq(i_neq_left,x(i,j2,f2))=1;
                    Bneq(i_neq_left,1)=1;
                    i_neq_left=i_neq_left+1;             
                end
            end
        end
    end
end
i_neq_right=i_neq_left;


%(16)
for i=1:1:n_path
    for f1=1:1:f_total
        for f2=1:1:f_total
            for j=1:1:l_max-1
                Aneq(i_neq_left,g1(f1,f2,j,j+1,i))=1;
                Aneq(i_neq_left+1,g1(f1,f2,j,j+1,i))=-1;
            end
            Aneq(i_neq_left,g2(f1,f2,i))=-l_max;
            Aneq(i_neq_left+1,g2(f1,f2,i))=1;
            i_neq_left=i_neq_left+2;
        end
    end
end
Bneq(i_neq_right:(i_neq_left-1),1)=0;
i_neq_right=i_neq_left;

%(20)
for i=1:1:n_path
    for f1=1:1:f_total
        for f2=1:1:f_total
            Aneq(i_neq_left,g2(f1,f2,i))=1;
            Aneq(i_neq_left,delta_ij(f1,f2))=-1;
            i_neq_left=i_neq_left+1;
        end
    end
end
Bneq(i_neq_right:(i_neq_left-1),1)=0;
i_neq_right=i_neq_left;

save('common_2_1.mat','-v7.3','Aneq','Bneq','i_neq_left','i_neq_right');

function index=x(i,j,f)
    global l_max
    global f_total
    index=(i-1)*l_max*f_total+(j-1)*f_total+f;
end

function index=delta_k(k)
    global n_path
    global l_max
    global f_total
    index=n_path*l_max*f_total+k;
end

function index=delta_ij(f1,f2)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+(f1-1)*f_total+f2;
end

function index=g1(f1,f2,j1,j2,i)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        (f1-1)*f_total*l_max*l_max*n_path+(f2-1)*l_max*l_max*n_path+ ...,
        (j1-1)*l_max*n_path+(j2-1)*n_path+i;
end

function index=g2(f1,f2,i)
    global n_path
    global l_max
    global f_total
    global n_dummy
    index=n_path*l_max*f_total+n_dummy+f_total*f_total+ ...,
        f_total*f_total*l_max*l_max*n_path+ ...,
        (f1-1)*f_total*n_path+(f2-1)*n_path+i;
end

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
    
    