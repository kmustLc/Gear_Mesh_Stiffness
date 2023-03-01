function h=symdisp(s)
%//SYMDISP Display a symbolic expression in human readable form.
%// symdisplay(S) displays the symbolic expression S in a small figure window,
%// using standard mathematical notation.
%// 
%// Examples:
%//   syms x t positive
%//   f=taylor(cos(x));
%//   symdisp(f)
%//   f=int(exp(-t)*t^(x-1),t,0,inf);
%//   symdisp(f)
%//
%// Required toolbox: Symbolic Math
%//
%// See also SYMBOLIC PRETTY.
if ~isa(s,'sym')
    s=sym(s);
    %error('�������������sym����,��ʹ�� sym() ����Ľ��ת��Ϊsym����.')
end
S=['$',latex(s),'$'];
S=strrep(S,'&','& \quad');
S=strrep(S,'{\it','\mathrm{');
h=msgbox(S,'�ַ�����ѧչʾ��ʽ');
h1=get(h,'children');
h2=h1(1);
h3=get(h2,'children');
if isempty(h3)
    h2=h1(2); h3=get(h2,'children');
end
set(h3,'visible','off')
set(h3,'interpreter','latex')
set(h3,'string',S)
set(h3,'fontsize',20)
w=get(h3,'extent');
W=get(h,'position');
W(3)=max(w(3)+10,125);
W(4)=w(4)+40;
set(h,'position',W)
h4=h1(2);
if ~strcmp(get(h4,'tag'),'OKButton'), h4=h1(1); end
o=get(h4,'position');
o(1)=(W(3)-o(3))/2;
set(h4,'position',o)
set(h3,'visible','on')
set(h,'color','w');
