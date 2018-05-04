function text_below(x,y,scale,corx,fontsize,ttt)
x=x-2;
text(scale*(x-corx),scale*(y+0.2+0.5),ttt,'HorizontalAlignment','center','VerticalAlignment','top','color','k','FontSize',fontsize);
% line([x-0.5 x x+0.5],y+[0 -0.2 0]+0.5,'color','k','LineWidth',0.5);
line(scale*[x-0.5 x x+0.5]-corx,scale*(y+[-0.0 +0.2 -0.0]+0.5),'color','k','LineWidth',2/4);

end