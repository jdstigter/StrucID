function xdot = xdotTemplateInputVar(t,x,u,th,f,vAlg)
v=vAlg(x,th,u);
xdot=f(t,x,u,th,v);