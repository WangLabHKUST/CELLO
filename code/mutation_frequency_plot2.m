function h = mutation_frequency_plot2(MutMat, GeneNames)
% conjoint distribution of two binomial distribution of proportion test
% Input matrix MutMat is an M * 4 matrix:
%       M total number of mutations
%       For each mutation, we have four numbers [x1 N1 x2 N2]:
%           x1: number of mutations in type1
%           N1: total number of cases in type1
%           x2: number of mutations in type2
%           N2: total number of cases in type2

scrsz = get(0,'Screensize');
h = figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)],'visible','on');
hold on

if nargin == 0
    disp('No input! Use CLL data from Hossein as an example.')
    MutMat=[16    58     18    174
        14    58     22    174
        11    58     12    174
        6    58      2    174    ];
    GeneNames = {'TP53','NOTCH1','SF3B1','FAT1'};
    M = 4;
elseif nargin == 1
    [ M, tmp_e4 ] = size( MutMat );
    if M > 14 || tmp_e4 ~= 4 || any( MutMat( : , 2 ) < MutMat( : , 1 ) ) || any( MutMat( : , 4 ) < MutMat( : , 3 ) )
        error( 'No Gene Name mode: Input error!' )
    else
        for i  = 1 : M
            GeneNames{ i } = ['gene' num2str( i )];
        end
    end
elseif nargin > 2
    error( 'Too many input arguments!' )
else
    [ M, tmp_e4 ] = size( MutMat );
    if M > 14 || tmp_e4 ~= 4 || any( MutMat( : , 2 ) < MutMat( : , 1 ) ) || any( MutMat( : , 4 ) < MutMat( : , 3 ) ) || M ~= numel( GeneNames )
        error( 'Gene Name mode: Input error!' )
    end
end

cutoff1 = 0.5;
cutoff2 = 0.65;
cutoff3 = 0.95;

N_circle = 1;

basiccolor='rgbmcpolasvfny'; % maximum 14 genes!!!
map=vivid(M*(N_circle+3),basiccolor(1:M),[0.5 1]); % light color based on rgbm

for m = 1 : M
    x1=MutMat(m,1);
    N1=MutMat(m,2);
    x2=MutMat(m,3);
    N2=MutMat(m,4);
    p1 = x1/N1;
    p2 = x2/N2;
    
    p1tmp = binopdf( 0:N1 , N1 , p1 ); % density of dimension 1
    p2tmp = binopdf( 0:N2 , N2 , p2 ); % density of dimension 2
    PMATRIX = p1tmp' * p2tmp; % joint density (assume two dimension are independent)
    
    all_p = sort(reshape(PMATRIX,(N1+1)*(N2+1),1),'descend');
    cut_density_p1 = all_p( find(cumsum(all_p)>=cutoff1, 1) );
    cut_density_p2 = all_p( find(cumsum(all_p)>=cutoff2, 1) );
    cut_density_p3 = all_p( find(cumsum(all_p)>=cutoff3, 1) );
    
    cut_density_p = [cut_density_p3 cut_density_p1 cut_density_p2];
    
    max_density_p = all_p(1);
    
    x=repmat( (0:N2)/N2 , N1+1 , 1 );
    y=repmat( ((0:N1)')/N1 , 1 , N2+1 );
    
    [~,h(m)] = contour(x,y,PMATRIX,'LineWidth',1,'LineStyle','--');
    
    hold on
    
     set( h(m) , 'ShowText','off' ,...
         'LevelList',cut_density_p(3):(max_density_p-cut_density_p(3))/N_circle:max_density_p)
     
     Cld = get(h(m), 'Children');
     for j=1:length(Cld)
         if strcmp(get(Cld(j), 'Type'), 'patch')
            set(Cld(j),'EdgeColor',map((m-1)*(N_circle+3)+j,:));
         end
     end
    
end

line('LineStyle','--','linewidth',2,'color',[0.5 0.5 0.5])

set( gca , 'xlim' , [-0.05,0.5] ,'ylim' , [-0.05 ,0.5] ,'fontsize',16 , 'xtick' , 0:0.1:1 , 'ytick',0:0.1:1)

for m = 1 : M
    x1=MutMat(m,1);
    N1=MutMat(m,2);
    x2=MutMat(m,3);
    N2=MutMat(m,4);
    plot( x2/N2 , x1/N1 , 'o' , 'color' , map((m-1)*(N_circle+3)+1,:) ,'linewidth',10 , 'MarkerSize' , 10 )
    text(x2/N2 , x1/N1,['',GeneNames{m}],'FontName','Arial','FontSize',18,'color',map((m-1)*(N_circle+3)+1,:))
end

set(gca,'LineWidth',2,'tickdir','out','FontSize',20,'FontWeight','bold')
axis square
