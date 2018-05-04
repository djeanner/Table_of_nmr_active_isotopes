% 0: Z/A
% 1: Z/N
% 2: N/Z
% add mass on the left

%
clear all
close all
scale=1;
fast=0;
shiftx=scale*60;% this is how much shift second part of n/z chart
%for switch_shift_y=[0 1 2]
% % for resol=[0 600]
% % for loop_unstable_isotopes=[ 0:1]
% %     for switch_shift_y=[ 0:2]
<<<<<<< HEAD
for loop_unstable_isotopes=[  1 0]%0  1
=======
for loop_unstable_isotopes=[0  1 ]%0  1
>>>>>>> 3ea4f53ea0351019a9b54a493db861e60768511c
    for switch_shift_y=[0 1 2]% 0 1 2
        clear pos_2d
        clear pos_2d_non_aliased
        disp('Initialize')
        if switch_shift_y==0
            txth='Z';
            txtv='A';
            
            type='Z_A_chart';
            fontsize=3.5;
            fontsize2=fontsize;
            
            maxypos=64;%nb cell vertical... used for modulo
            
            tab_shiftx=[0 10+2 20 30-1 0];
           maxypos=51;%nb cell vertical... used for modulo
           tab_shiftx=6*[0 1 2 3 4 5 6 7];
           tab_shiftx=[0 5+1 10 15-1 20-5 25-8 0];
            
        end
        if switch_shift_y==1
            txth='Z';
            txtv='N';
            type='Z_N_chart';
            
            fontsize=3.5;
            
            fontsize2=fontsize-0.5;
            
            
            maxypos=66;%nb cell vertical... used for modulo
            tab_shiftx=[0 16 32 0 ];
            maxypos=42;%nb cell vertical... used for modulo
            tab_shiftx=[0 6 12-6 18-10 0];
             maxypos=48;%nb cell vertical... used for modulo
            tab_shiftx=[0 8+1 16+1 24+1 32+1 0];
            
        end
        if switch_shift_y==2
              txth='N';
            txtv='Z';
            type='N_Z_chart';
            
            fontsize=1.5;
            fontsize=3.5;
            
            fontsize2=fontsize-0.3;
            fontsize2=fontsize-0.5;
            
            maxypos=64;%nb cell vertical... used for modulo
                                                            tab_shiftx=[0 60 0 ];

        end

        figure(switch_shift_y+1);clf;drawnow;hold on
        tic
        if switch_shift_y==0
                    origin=[-0.5 0.5];
        else
                    origin=[-0.5 -0.5];
        end
                    length_arrows_origin=4;
                    oi=0.5/2;
                    oa=0.2;
                    lwar=1;
                    %line of arrow
                    plot([origin(1,1) origin(1,1)],[origin(1,2) origin(1,2)+length_arrows_origin],'k-','LineWidth',lwar*scale);hold on
                    plot([origin(1,1) origin(1,1)+length_arrows_origin ],[origin(1,2) origin(1,2)],'k-','LineWidth',lwar*scale);hold on
                   %tip arrow
                    plot([origin(1,1)+oi origin(1,1) origin(1,1)-oi],length_arrows_origin+origin(1,2)+[ -oa 0 -oa],'k-','LineWidth',lwar*scale);hold on
                    plot(length_arrows_origin+origin(1,2)+[ -oa 0 -oa],[origin(1,1)+oi origin(1,1) origin(1,1)-oi],'k-','LineWidth',lwar*scale);hold on
%text

little_space_between_top_arrow_and_txt=0.1;
fontsizea=fontsize*1.2;
    text(origin(1,1)+length_arrows_origin+little_space_between_top_arrow_and_txt,origin(1,2),txth,'HorizontalAlignment','left','VerticalAlignment','middle','color','k','FontSize',fontsizea)
if switch_shift_y==2
    text(origin(1,1),origin(1,2)+length_arrows_origin+little_space_between_top_arrow_and_txt,txtv,'HorizontalAlignment','center','VerticalAlignment','bottom','color','k','FontSize',fontsizea)
else
      text(origin(1,1),origin(1,2)+length_arrows_origin+little_space_between_top_arrow_and_txt,txtv,'HorizontalAlignment','center','VerticalAlignment','top','color','k','FontSize',fontsizea)
  
end
        set(gcf,'color','w');
        disp('Read data')
        toc
        s=tdfread('./NMR_data.txt');
        Per_tab={'n','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br',...
            'Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',...
            'Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',...
            'Rf'};
        if switch_shift_y<2
            pos_x_table_above=10000*ones(4,size(Per_tab,2));
            pos_x_table_below=-10000*ones(4,size(Per_tab,2));
        else
            pos_x_table_above=10000*ones(4,size(s.A,1));
            pos_x_table_below=-10000*ones(4,size(s.A,1));
            
        end
        %if switch_shift_y==0
        pos_iso=10000*ones(1,1+max(max(s.N-s.A*00000000)));
        pos_iso2=zeros(1,1+max(max(s.N-s.A*00000000)));
        %else
        %    pos_iso=10000*ones(1,1+max(max(s.N-s.A*switch_shift_y)));
        
        %end
        for loopi=1:size(s.A,1)
            txt=s.spin0x28halfs0x29(loopi,:);
            if strcmp(txt,'#VALUE!')
                s.Nu(loopi)=0;
                spintab(loopi)=1e6;
            else
                tmp_value=str2num(txt);
                if size(tmp_value,2)>0
                    spintab(loopi)=tmp_value;
                else
                    spintab(loopi)=1e6;
                    
                end
            end
        end
        disp('test data')
        toc
        
        
        s.Nu(585)=11.4;%???
        lim=2.7;
        %lim=100002.7;
        for loopi=1:size(s.Nu,1)
            if abs(s.Nu(loopi))>lim
                s.Nu(loopi)=lim*sign(s.Nu(loopi));
                
            end
        end
        if fast
            nb_colors=16;
        else
            nb_colors=512;
            
        end
        com_col=colormap(flipud(summer(nb_colors)));
        % com_col=colormap(copper(64));
        %  com_col=colormap('hsv');
        for loopi=1:size(com_col,1)
            a1=maxypos*(size(com_col,1)-loopi)/size(com_col,1);
            a2=maxypos*(size(com_col,1)-loopi+1)/size(com_col,1);
            if switch_shift_y~=2
                filcol=com_col(loopi,:);
            else
                filcol=com_col(size(com_col,1)-loopi+1,:);
                
            end
            fill(scale*([-1 -1 0 0 -1]-1.5*0+105-3),scale*[a1 a2 a2 a1 a1],filcol,'EdgeColor','none','LineWidth',1*scale);
            
            
            
        end
        drawnow
        
        disp('test data2')
        toc
        %search max(abs(nu))
        tmpdelll=s.Nu';
        maxnu=max(max(abs(tmpdelll./spintab)));
        table_AN=zeros(max(max(s.A))+1,max(max(s.N))+1);
        ignore=zeros(1,size(s.A,1));
        
        %make table of spins to ignore (keep the most stable when
        %multiple)
        disp(['Loop over ' num2str(size(s.A,1)) ' isotopes/spin'])
        
        for loopi=1:size(s.A,1)
            
            % check multiplet spins for given isotopes...
            tA=s.A(loopi,1)+1;
            tN=s.N(loopi,1)+1;
            loopi2=loopi;
            if  table_AN(tA,tN)~=0
                
                disp(['found multiple spins for isotope ' num2str(tN-1) '' Per_tab{tA} ' Keep longer lived'])
                
                
                %set t 1/2
                tmptxt=s.T10x2F2(loopi,:);
                del='#VALUE!    ';
                test=strcmp(tmptxt,del);
                if test
                    t1_2=-1;
                else
                    if strcmp(tmptxt,'')
                        t1_2=-2;
                    else
                        t1_2=str2num(tmptxt);
                    end
                end
                cur_t1_2=t1_2;
                %    disp(['Current (' num2str(loopi) '): t1/2 = ' num2str(cur_t1_2)])
                %set t 1/2
                tmptxt=s.T10x2F2(table_AN(tA,tN),:);
                del='#VALUE!    ';
                test=strcmp(tmptxt,del);
                if test
                    t1_2=-1;
                else
                    if strcmp(tmptxt,'')
                        t1_2=-2;
                    else
                        t1_2=str2num(tmptxt);
                    end
                end
                if cur_t1_2<=0 cur_t1_2=1e50;end
                if t1_2<=0 t1_2=1e50;end
                %     disp(['other (' num2str(table_AN(tA,tN)) '): t1/2 = ' num2str(t1_2)])
                if cur_t1_2<t1_2
                    ignore(1,loopi)=1;
                    loopi2=table_AN(tA,tN);
                    
                else
                    ignore(1,table_AN(tA,tN))=1;
                end
            end
            table_AN(tA,tN)=loopi2;
            
        end
        toc
        disp('second step')
        disp(['Loop over ' num2str(size(s.A,1)) ' isotopes/spin'])
        for loopi=1:size(s.A,1)
            if loopi>50
                %   drawnow
                %   break
            end
            if ignore(1,loopi)>0
                tA=s.A(loopi,1)+1;
                tN=s.N(loopi,1)+1;
                disp(['Ignored a spin for ' num2str(tN-1) '' Per_tab{tA} ' (' num2str(loopi) ')' ])
            else
                %set parameters
                tmptxt=s.T10x2F2(loopi,:);
                del='#VALUE!    ';
                test=strcmp(tmptxt,del);
                if test
                    t1_2=-1;
                else
                    if strcmp(tmptxt,'')
                        t1_2=-2;
                    else
                        t1_2=str2num(tmptxt);
                    end
                end
                %set nu
                nu=s.Nu(loopi);
                %set spin text
                txt=s.spin0x28halfs0x29(loopi,:);
                if strcmp(txt,'#VALUE!')
                    txt='';spin=[];
                else
                    spin=str2num(txt);
                    if spin/2==round(spin/2)
                        txt=num2str(spin/2);
                    else
                        txt=[num2str(spin) '/2'];
                    end
                end
                if switch_shift_y<2
                    
                    x=s.A(loopi,1);
                    
                    y=s.N(loopi,1);
                else
                    x=s.N(loopi,1);
                    
                    y=s.A(loopi,1);
                end
                y_orig=y;
                x_orig=x;
                if switch_shift_y==1
                    y=y-x;
                end
                corx=0;
                corxn=0;
                if switch_shift_y==2
                    x=x-y;
                    yo=y;
                    y=mod(y,maxypos);
                                                            shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));

                    if y==yo
                        
                        corx=0;
                    else
                        corx=scale*shiftx;
                        
                    end
                else
                    yo=y;
                    y=mod(y,maxypos);
                    
                    shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));
                    corx=scale*shiftx;
                     if (-1+ 1+round((yo-y)/maxypos))>0
                      shiftx=tab_shiftx(1, -1+ 1+round((yo-y)/maxypos));
                     end
                    corxn=scale*shiftx;
                    
                    
                end
                if size(spin,2)>0
                    pos1_64=1+round((nb_colors-1)*abs(nu/spin)/maxnu);
                    % if pos1_64>40
                    % pos1_64
                    % end
                    filcol1='g';
                    
                    filcol1=com_col(pos1_64,:);
                    
                    if (t1_2==0) || (t1_2>(100*365*24*3600)) % larger than 100 years
                        col='k';
                        filcol= filcol1;
                        coltxt='w';
                        coltxt='k';
                        
                    else
                        col='k';
                        col=0.8*[1 1 1];
                        filcol= filcol1;
                        coltxt='k';
                        
                        filcol= filcol1*0.3+0.7;
                        filcol= 'w';
                        
                    end
                    unstable_isotope=(t1_2<1e-3) && (t1_2>0);
                    if ~(unstable_isotope && loop_unstable_isotopes) % this is to eliminate the nuclear excited states
                        line_width=scale*0.5;
                        
                        color_frame='k';
                        disx=0.5;disy=0.5;
                        
                        
                        fill(scale*[x+disx x+disx x-disx x-disx x+disx]-corx,scale*[y+disy y-disy y-disy y+disy y+disy],filcol,'LineWidth',line_width)
                        if line_width>0
                            line(scale*[x+disx x+disx x-disx x-disx x+disx]-corx,scale*[y+disy y-disy y-disy y+disy y+disy],'LineWidth',line_width,'color',color_frame')
                        end
                        if unstable_isotope
                            color_frame='r';
                            disx=0.45;disy=disx;
                            line_width=scale*0.3;
                            if line_width>0
                                line(scale*[x+disx x+disx x-disx x-disx x+disx]-corx,scale*[y+disy y-disy y-disy y+disy y+disy],'LineWidth',line_width,'color',color_frame')
                            end
                        end
                        text(scale*x-corx,scale*y,txt,'HorizontalAlignment','center','FontSize',fontsize,'color',coltxt);
                        
                        
                        if pos_x_table_above(1,x+1)>y
                            if pos_x_table_above(1,x+1)>1000
                                pos_x_table_above(1,x+1)=y;
                                                        pos_x_table_above(3,x+1)=corx;
                                                        
                            else
                                if abs(pos_x_table_above(1,x+1)-y)>(maxypos/2) % store two values...
                                    pos_x_table_above(2,x+1)=pos_x_table_above(1,x+1);
                                    pos_x_table_above(4,x+1)=corxn;
                                    
                                end
                                pos_x_table_above(1,x+1)=y;
                                pos_x_table_above(3,x+1)=corx;
                                
                            end
                        end
                        if pos_x_table_below(1,x+1)<1000
                            pos_x_table_below(1,x+1)=y;
                            pos_x_table_below(3,x+1)=corx;
                        end
                        
                        if abs(pos_x_table_below(1,x+1)-y)>(maxypos/2) % store two values...
                            if pos_x_table_below(2,x+1)<y
                                pos_x_table_below(2,x+1)=y;
                                pos_x_table_below(4,x+1)=corx;
                            end
                        else
                            if pos_x_table_below(1,x+1)<y
                                pos_x_table_below(1,x+1)=y;
                                pos_x_table_below(3,x+1)=corx;
                            end
                            
                            %   pos_x_table_below(2,x+1)=pos_x_tapos_x_table_belowble_above(1,x+1);
                            %  pos_x_table_below(4,x+1)=corxn;
                            
                        end
%                                 pos_x_table_below(1,x+1)=y;
%                                 pos_x_table_below(3,x+1)=corx;
%                                 
%                             end
%                         end
                        if pos_iso(1,y_orig+1)>x
                            pos_iso(y_orig+1)=x;
                        end
                        if pos_iso2(1,y_orig+1)<x
                            pos_iso2(y_orig+1)=x;
                        end
                        %                         if y==22
                        %                             [x y]
                        %                             plot(x,y,'ro')
                        %                         end
                        pos_2d(x+2,y+2)=1;
                        if switch_shift_y==2
                            
                            pos_2d_non_aliased(x+2,yo+2)=1;
                        end
                        
                    end
                end
                if (loopi==1)
                    axis off
                    
                    %update...
                    
                    axis(scale*[-1.5 105 -1 maxypos+1])
                    
                    if switch_shift_y<2
                        
                        set(gca,'Ydir','reverse')
                    end
                    axis('equal')
                end
                if mod(loopi,100)==0
                    
                    drawnow
                end
            end
        end
        drawnow
        
        disp('Display element')
        toc
        %display Element
        if switch_shift_y<2
            
            for looopi=1:size( pos_x_table_above,2)
                x=looopi-1;
                for indi=1:2
                 y=(pos_x_table_above(indi,looopi)-1)';
                    corx=pos_x_table_above(2+indi,looopi);
%                                                             shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));
% 
%                 if y==yo
%                     corx=0;
%                 else
%                     corx=scale*shiftx;
%                 end
                  %   corx=0;

                  if y<1000
                      text_above(x+2,y+1,scale,corx,fontsize,Per_tab{looopi})
                      dont_write(x+2,y+2)=1;% this is to avoid writing isotop text on the element text
                      pos_2d(x+2,y+2)=2;
                      
                      
                  end
                end
                for indi=1:2
                 y=(pos_x_table_below(indi,looopi)-1)';
                    corx=pos_x_table_below(2+indi,looopi);
%                                                             shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));
% 
%                 if y==yo
%                     corx=0;
%                 else
%                     corx=scale*shiftx;
%                 end
                  %   corx=0;

                  if y>-1000
                    
                      
                      
                      text_below(x+2,y+1,scale,corx,fontsize,num2str(looopi-1))
                    %  if y>0
                      dont_write(x+2,y+4*1)=1;% this is to avoid writing isotop text on the element text
                      pos_2d(x+2,y+4*1)=2;
                   %   end
                      
                  end
                end
            end
        end
        if switch_shift_y==2
            
            for looopi=1:size( Per_tab,2)
                
                yo=looopi-1;
                y=(mod(looopi-1,maxypos));
                                                            shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));

                if y==yo
                    corx=0;
                else
                    corx=scale*shiftx;
                end
                x=pos_iso(looopi)-1;
                if x<1000
                    % text(x,y-0.15+0.45,Per_tab{looopi},'HorizontalAlignment','center','VerticalAlignment','bottom','color','k','FontSize',fontsize);
                    text(scale*(x-0.2*0+0.15)-corx,scale*y,Per_tab{looopi},'HorizontalAlignment','right','VerticalAlignment','middle','color','k','FontSize',fontsize);
                    %  line([x-0.5 x x+0.5],y+[-0.0 -0.2 -0.0]+0.5,'color','k','LineWidth',0.5);
                    line(scale*(x+[-0.0 -0.2 -0.0]+0.5)-corx,scale*[y-0.5 y y+0.5],'color','k','LineWidth',0.5);
                    dont_write(x+2,y+2)=1;% this is to avoid writing isotop text on the element text
                    pos_2d(x+2,y+2)=2;
                    pos_2d_non_aliased(x+2,y+2)=2;
                    
                end
            end
        end
        drawnow
        disp('display isotope')
        
        toc
        dont_write(size(dont_write,1)+5,size(dont_write,2)+15)=0;%increase size or array to avoid zeros...
        pos_2d(size(pos_2d,1)+5,size(pos_2d,2)+15)=0;%increase size or array to avoid zeros...
        for looopi=1:size( pos_iso,2)
            
            x=pos_iso(looopi);
            
            if x<1000
                if switch_shift_y==0
                    yo=looopi-1;
                    y=1+mod(looopi-1,maxypos);
                      shiftx=tab_shiftx(1,1+round((yo+1-y)/maxypos));

                if y==yo
                    corx=0;
                else
                    corx=scale*shiftx;
                end
                    if pos_2d(x+2-1,y+2-1)==0
                        
                        text(scale*(x-0.5-0.25)-corx,scale*(y-1),num2str(looopi-1),'HorizontalAlignment','right','VerticalAlignment','middle','color','k','FontSize',fontsize2);
                        line(scale*([x x-0.2 x]-0.5)-corx,scale*(y+[-0.5 -0.0 0.5]-1),'color','k','LineWidth',0.5);
                    else
                        disp(['no isotop label for ' num2str(looopi-1) 'because already somehting in this place:' num2str(pos_2d(x+2-1,y+2-1)) ' ' num2str(x) ' ' num2str(y)])
                    end
                end
                if switch_shift_y==1
                    
                    y=1+mod(looopi-1-x,maxypos);
                    yo=looopi-1-x;
                    shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));

                if y==yo
                    corx=0;
                else
                    corx=scale*shiftx;
                end
                    %bottom left
                    if dont_write(x+2 -1,y+2 )==0
                        
                        text(scale*(x-0.5-0.25+0.1)-corx,scale*(y-1+1.15),num2str(looopi-1),'HorizontalAlignment','right','VerticalAlignment','bottom','color','k','FontSize',fontsize2);
                        line(scale*([x x-0.2 ]-0.5)-corx,scale*(y+[0.5 0.5+0.2]-1),'color','k','LineWidth',0.5);
                    end
                    %top right
                    x=pos_iso2(looopi);
                    y=1+mod(looopi-1-x,maxypos);
                     yo=looopi-1-x;
                    shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));

                if y==yo
                    corx=0;
                else
                    corx=scale*shiftx;
                end
                    if dont_write(x+2 +1,y+2 -2)==0
                        
                        text(scale*(x+0.5+0.25-0.1)-corx,scale*(y-1-1.15),num2str(looopi-1),'HorizontalAlignment','left','VerticalAlignment','top','color','k','FontSize',fontsize2);
                        line(scale*([x x+0.2 ]+0.5)-corx,scale*(y+[-0.5 -0.5-0.2]-1),'color','k','LineWidth',0.5);
                    end
                end
                
            end
        end
        if switch_shift_y==2
            
            for looopi=1:size( pos_2d_non_aliased,1)-2
                
                x=(looopi);
                tmp=pos_2d_non_aliased(x+2,:);
                tt=find(tmp==1);
                if size(tt,2)>0
                    %%% above
                    pos=size(tt,2);
                    yo=tt(1,pos)-2;
                    
                    y=mod(yo-0,maxypos)+0;
                                                            shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));
                    if y==yo
                        corx=0;
                    else
                        corx=scale*shiftx;
                    end
                    if pos_2d_non_aliased(x+2,y+2+1)==0
                        text(scale*(x-corx),scale*(y   +0.15+0.45),num2str(looopi-1+1),'HorizontalAlignment','center','VerticalAlignment','bottom','color','k','FontSize',fontsize);
                        % line([x-0.5 x x+0.5],y+[0 -0.2 0]+0.5,'color','k','LineWidth',0.5);
                        line(scale*[x-0.5 x x+0.5]-corx,scale*(y    -[-0.0 -0.2 -0.0]+0.5),'color','k','LineWidth',0.5);
                        
                        %      text(scale*(x-0.5-0.25+0.1)-corx,scale*(y-1+1.15),num2str(looopi-1),'HorizontalAlignment','right','VerticalAlignment','bottom','color','k','FontSize',fontsize2);
                        %     line(scale*([x x-0.2 ]-0.5)-corx,scale*(y+[0.5 0.5+0.2]-1),'color','k','LineWidth',0.5);
                    end
                    %%% below
                    pos=1;
                    yo=tt(1,pos)-2;
                    
                    y=mod(yo-1,maxypos)+1;
                    if yo>0
                                                            shiftx=tab_shiftx(1,1+round((yo-y)/maxypos));
                    else
                       shiftx=0; 
                    end

                    if y==yo
                        corx=0;
                    else
                        corx=scale*shiftx;
                    end
                    if pos_2d_non_aliased(x+2,y+2-1)==0
                        text(scale*(x-corx),scale*(y   -0.15-0.45),num2str(looopi-1+1),'HorizontalAlignment','center','VerticalAlignment','top','color','k','FontSize',fontsize);
                        % line([x-0.5 x x+0.5],y+[0 -0.2 0]+0.5,'color','k','LineWidth',0.5);
                        line(scale*[x-0.5 x x+0.5]-corx,scale*(y    +[-0.0 -0.2 -0.0]-0.5),'color','k','LineWidth',0.5);
                        
                        %      text(scale*(x-0.5-0.25+0.1)-corx,scale*(y-1+1.15),num2str(looopi-1),'HorizontalAlignment','right','VerticalAlignment','bottom','color','k','FontSize',fontsize2);
                        %     line(scale*([x x-0.2 ]-0.5)-corx,scale*(y+[0.5 0.5+0.2]-1),'color','k','LineWidth',0.5);
                    end
                end
            end
        end
        
        drawnow
        disp('Rescale')
        toc
        drawnow
        
        
            axis(scale*[-1.5 105 -1 maxypos+1])
          %  axis(scale*[-5 10 -5 10])
            
        
        axis off
        if switch_shift_y<2
            
            set(gca,'Ydir','reverse')
        end
        axis('equal')
        
        si=[ 42.0 29.7 ]-1*2;%for A3 with 1 cm margin at both sides
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf, 'PaperSize', si);
        set(gcf, 'PaperPositionMode', 'manual');
        set(gcf, 'PaperPosition', [1 1 si]);
        fig = gcf;
        %  fig.InvertHardcopy = 'off';%to avoid change color
        %set(findall(gcf,'-property','FontSize'))
        %set(findall(gcf,'-property','FontName'),'FontName','Times New Roman')
        
        % Export file
        
        %  print('-r1200','-dpdf',['test' num2str(switch_shift_y) '_1200dpi.pdf'])
        if ~exist(['./' type 's'],'dir')
            mkdir(['./' type 's']);
        end
        disp('dump file')
        
        toc
        drawnow
        for resol=[600-fast*400  ]% 600 %0 for svg 300 (fast) 600 for high resolution pdf
            
            if loop_unstable_isotopes
                % print('-r600','-dpdf',['test' num2str(switch_shift_y) '_600dpi.pdf'])
                if resol==0
                    plot2svg(['./' type 's/' type  '.svg'])
                else
                    print(['-r' num2str(resol)],'-dpdf',['./' type 's/' type '_' num2str(resol) 'dpi.pdf'])
                end
            else
                if resol==0
                    plot2svg(['./' type 's/' type '_with_unstable_isotopes.svg'])
                else
                    print(['-r' num2str(resol)],'-dpdf',['./' type 's/' type '_with_unstable_isotopes_' num2str(resol) 'dpi.pdf'])
                end
            end
            %   print('-r600','-dtiff',['test' num2str(switch_shift_y) '_600dpi.tif'])
        end
        disp('end of dump file')
        toc
          %      fzuk

    end
end

