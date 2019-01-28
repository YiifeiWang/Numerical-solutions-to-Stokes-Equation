function [U,V,P] = DGS_InnerM(N, U0, V0, P0)
	U = U0; V = V0; P = P0;
	u = reshape(U,N-1,N);
	v = reshape(V,N-1,N)';
	p = reshape(P,N,N);
	h = 1/N;
	for i=2:N-1
		for j=2:N-1
			r = (u(i,j)-u(i-1,j)+v(i,j)-v(i,j-1))/h;
			delta = r*h/4;
			u(i-1,j) = u(i-1,j)+delta;
			u(i,j) = u(i,j)-delta;
			v(i,j-1) = v(i,j-1)+delta;
			v(i,j) = v(i,j)-delta;
			p(i,j) = p(i,j)-r;
			p(i+1,j) = p(i+1,j)+r/4;
			p(i-1,j) = p(i-1,j)+r/4;
			p(i,j+1) = p(i,j+1)+r/4;
			p(i,j-1) = p(i,j-1)+r/4;
		end
	end

	i = 1;
    for j=2:N-1
        r = (u(i,j)+v(i,j)-v(i,j-1))/h;
        delta = r*h/3;
        u(i,j) = u(i,j)-delta;
        v(i,j-1) = v(i,j-1)+delta;
        v(i,j) = v(i,j)-delta;
        p(i,j) = p(i,j)-4*r/3;
        p(i+1,j) = p(i+1,j)+r/3;
        p(i,j+1) = p(i,j+1)+r/3;
        p(i,j-1) = p(i,j-1)+r/3;
    end

    i = N;
    for j=2:N-1
        r = (-u(i-1,j)+v(i,j)-v(i,j-1))/h;
        delta = r*h/3;
        u(i-1,j) = u(i-1,j)+delta;
        v(i,j-1) = v(i,j-1)+delta;
        v(i,j) = v(i,j)-delta;
        p(i,j) = p(i,j)-4*r/3;
        p(i-1,j) = p(i-1,j)+r/3;
        p(i,j+1) = p(i,j+1)+r/3;
        p(i,j-1) = p(i,j-1)+r/3;
    end

    j = 1;
    for i=2:N-1
        r = (u(i,j)-u(i-1,j)+v(i,j))/h;
        delta = r*h/3;
        u(i-1,j) = u(i-1,j)+delta;
        u(i,j) = u(i,j)-delta;
        v(i,j) = v(i,j)-delta;
        p(i,j) = p(i,j)-4*r/3;
        p(i+1,j) = p(i+1,j)+r/3;
        p(i-1,j) = p(i-1,j)+r/3;
        p(i,j+1) = p(i,j+1)+r/3;
    end

    j = N;
    for i=2:N-1
        r = (u(i,j)-u(i-1,j)-v(i,j-1))/h;
        delta = r*h/3;
        u(i-1,j) = u(i-1,j)+delta;
        u(i,j) = u(i,j)-delta;
        v(i,j-1) = v(i,j-1)+delta;
        p(i,j) = p(i,j)-4*r/3;
        p(i+1,j) = p(i+1,j)+r/3;
        p(i-1,j) = p(i-1,j)+r/3;
        p(i,j-1) = p(i,j-1)+r/3;
    end

    i = 1; j = 1;
    r = (u(i,j)+v(i,j))/h;
    delta = r*h/2;
    u(i,j) = u(i,j)-delta;
    v(i,j) = v(i,j)-delta;
    p(i,j) = p(i,j)-2*r;
    p(i+1,j) = p(i+1,j)+r/2;
    p(i,j+1) = p(i,j+1)+r/2;

    i = 1; j = N;
    r = (u(i,j)-v(i,j-1))/h;
    delta = r*h/2;
    u(i,j) = u(i,j)-delta;
    v(i,j-1) = v(i,j-1)+delta;
    p(i,j) = p(i,j)-2*r;
    p(i+1,j) = p(i+1,j)+r/2;
    p(i,j-1) = p(i,j-1)+r/2;

    i = N; j = 1;
    r = (-u(i-1,j)+v(i,j))/h;
    delta = r*h/2;
    u(i-1,j) = u(i-1,j)+delta;
    v(i,j) = v(i,j)-delta;
    p(i,j) = p(i,j)-2*r;
    p(i-1,j) = p(i-1,j)+r/2;
    p(i,j+1) = p(i,j+1)+r/2;

    i = N; j = N;
    r = (-u(i-1,j)-v(i,j-1))/h;
    delta = r*h/2;
    u(i-1,j) = u(i-1,j)+delta;
    v(i,j-1) = v(i,j-1)+delta;
    p(i,j) = p(i,j)-2*r;
    p(i-1,j) = p(i-1,j)+r/2;
    p(i,j-1) = p(i,j-1)+r/2;


	U = u(:);
	vt = v';
	V = vt(:);
	P = p(:);

end