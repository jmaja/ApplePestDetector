function h=display_network(A,S_var)

[L M]=size(A);


figure(1)

dim=ceil(sqrt(M));
sz=sqrt(L);

h=zeros(M,1);

for i=1:M
  a_max=max(abs(A(:,i)));
  subplot(dim,dim,i)
  h(i)=imagesc(reshape(A(:,i),sz,sz),'EraseMode','none',[-a_max a_max]);
  axis square, axis off
  drawnow
end


if exist('S_var','var')
  figure(2)
  subplot(211), bar(S_var), title('s variance')
  subplot(212), bar(sqrt(sum(A.*A))), title('basis norm (L2)')
end

drawnow
