function [ sol ] = polysplinefit( x,y,m,cond )
    %   This function fits a polynomial splines of order m to a given data (x,y).
    %   it is valid for one dimension only.
    %   the function uses interpolation approach so it is not suitable for
    %   noisy data
    %   ----------------------------------------------
    %   inputs
    %   x           x data must be increasing. ex x = [1,2,3]
    %   y           f(x)
    %   m           polynomial order for each spline
    %   cond        the addition conditions needed to compute the splines
    %               it's given by spline initial derivative values
    %               for the first and last spline only.
    %               cond rows number = spline order - 1.
    %
    %   cond form is [spline,derivative_order,value_of_the_derivative]
    %                 spline : 1 for the first spline
    %                          2 for the last splines
    %   ex:
    %   cond = [1 1 0;2 1 0]
    %   this means that the value of the first derivative of the first spline
    %   is zero and the value of the first derivative of the last spline is zero
    %   
    %   ---------------------------------------------------------------------
    % 	output
    %   sol      	spline between each x
    %   ---------------------------------------------------------------------
    %   example
    %  x = [0 2 4 6 8];
    %  y = [0 2.2484 2.3164 2.5413 2.8626];
    %  m = 2;
    %  cond = [1 1 0;2 1 0;1 2 0;2 2 0];
    %   sol = polysplinefit(x,y,m,cond)
    %
    %   All copyrights goes to Mohammad Al-Fetyani
    %   University of Jordan
    if size(x,2) > size(x,1)
        x = x';
    end
    if size(y,2) > size(y,1)
        y = y';
    end
    nData = size(x,1); % number of data
    n = nData - 1;
    % Number of unknowns initialize
    nCoeff = (m+1)*n;
    a=zeros(nCoeff);
    b=zeros(nCoeff,1);
    %% boundary
    for j=0:m
        a(1,j+1)= x(1)^j;
        a(2*n,nCoeff-m+j)= x(nData)^j;
        b(1)=y(1);
        b(2*n)=y(nData);
    end
    %% position condition
    r=2;
    for i=2:n
        for j=0:m
            c = (m+1)*(i-1)+1;
            a(r,c-(m+1)+j)= x(i)^j;
            a(r+1,c+j)= x(i)^j;
        end
        b(r:r+1)=y(i);
        r=r+2;
    end
    %% derivative continuity conditions
    r=r+1;
    f = ones(1,m+1);
    df = 0;
    for i=m-1:-1:1
        f=polyder(f);
        df=df+1;
        co = 1;
        c = df+1;
        for j=2:n
            f1=flip(f).*x(co+1).^(0:i);
            a(r,c:i+c)=f1;
            d = c + i + 1 + df;
            a(r,d:i+d)=-f1;
            c =d;
            r=r+1;
            co=co+1;
        end
    end
    %% addtional constrain
    for i=r:nCoeff
        f=ones(1,m+1);
        j = 1+i-r;
        a1=zeros(1,nCoeff);
        u=m-cond(j,2);
        shift=0;
        for k=1:cond(j,2)
            f=polyder(f);
            shift=shift+1;
        end
        if cond(j,1) == 1
            f1=flip(f).*x(1).^(0:u);
            c=shift+1;
        else
            f1=flip(f).*x(end).^(0:u);
            c=nCoeff-(m)+shift;
        end
        a1(c:c+m-shift) = f1;
        a(i,:)=a1;
        b(r)=cond(j,3);
    end
    %% solution
    s=a\b;
    j=1;
    for i=1:m+1:nCoeff
        fn(j,:)={flip(s(i:i+m))',['[',num2str(x(j)),',',num2str(x(j+1)),']']};
        j=j+1;
    end
    sol = fn;
end
