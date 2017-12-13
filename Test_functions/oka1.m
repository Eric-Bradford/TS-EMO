function f = oka1( x )

	x1p = cos(pi./12.0).*x(:,1) - sin(pi./12.0).*x(:,2);
	x2p = sin(pi./12.0).*x(:,1) + cos(pi./12.0).*x(:,2);

	f(:,1) = x1p;
	f(:,2) = sqrt(2*pi) - sqrt(abs(x1p)) + 2 .* (abs(x2p-3.*cos(x1p)-3).^0.33333333);
end