function hcom = corrMut(kdstr)
% plot mut-mut correlations: co-occurrence vs. mutually exclusive in
% primary and recurrent tumor samples, respectively.

genelist = kdstr.kdlist;
ng = numel(genelist);

Pmat = full(spones(kdstr.Pmat + kdstr.Cmat))'; % mutation table of primary samples
Rmat = full(spones(kdstr.Rmat + kdstr.Cmat))'; % mutation table of recurrence samples

%%

tumortype = cell(ng*ng - ng, 1);
genex = cell(ng*ng - ng, 1);
geney = cell(ng*ng - ng, 1);
gi = zeros(ng*ng - ng, 1);
gj = zeros(ng*ng - ng, 1);
dotx = zeros(ng*ng - ng, 1); % scatter x-axis
doty = zeros(ng*ng - ng, 1); % scatter y-axis
sz = zeros(ng*ng - ng, 1); % scatter sizes
cr = zeros(ng*ng - ng, 3); % scatter colors
pval = ones(ng*ng - ng, 1);
odds = zeros(ng*ng - ng, 1);

pvcut = 0.1; % specified in ng3590

t = 0;
for i = 1:ng
    for j = 1:ng
        Z = zeros(2,2);
        if i < j % primary
            t = t + 1;
            tumortype{t} = 'Primary';
            genex{t} = genelist{i};
            geney{t} = genelist{j};
            gi(t) = i;
            gj(t) = j;
            dotx(t) = -i;
            doty(t) = ng - j + 1;
            Z(1,1) = nnz(Pmat(i,:) == 1 & Pmat(j,:) == 1);
            Z(1,2) = nnz(Pmat(i,:) == 1 & Pmat(j,:) == 0);
            Z(2,1) = nnz(Pmat(i,:) == 0 & Pmat(j,:) == 1);
            Z(2,2) = nnz(Pmat(i,:) == 0 & Pmat(j,:) == 0);
            [~,p,~] = fishertest(Z);
            pval(t) = p;
            % odds ratio
            Z1 = Z + ones(2,2);
            oz = Z1(1,1)*Z1(2,2)/(Z1(1,2)*Z1(2,1));
            odds(t) = oz;
            if p < pvcut
                sz(t) = -log10(p);
                if oz > 1
                    cr(t,:) = [1 0 0]; % red: co-mut
                else
                    cr(t,:) = [0 1 0]; % green: mutual exclusive
                end
            else
                sz(t) = 1;
                cr(t,:) = [0.75 0.75 0.75]; % grey
            end

        elseif i > j % relapsed
            t = t + 1;
            tumortype{t} = 'Recurrence';
            genex{t} = genelist{i};
            geney{t} = genelist{j};
            gi(t) = i;
            gj(t) = j;
            dotx(t) = j;
            doty(t) = ng - i + 1;
            Z(1,1) = nnz(Rmat(i,:) == 1 & Rmat(j,:) == 1);
            Z(1,2) = nnz(Rmat(i,:) == 1 & Rmat(j,:) == 0);
            Z(2,1) = nnz(Rmat(i,:) == 0 & Rmat(j,:) == 1);
            Z(2,2) = nnz(Rmat(i,:) == 0 & Rmat(j,:) == 0);
            [~,p,~] = fishertest(Z);
            pval(t) = p;
            % odds ratio
            Z1 = Z + ones(2,2);
            oz = Z1(1,1)*Z1(2,2)/(Z1(1,2)*Z1(2,1));
            odds(t) = oz;
            if p < pvcut
                sz(t) = -log10(p);
                if oz > 1
                    cr(t,:) = [0 0 0]; % black
                else
                    cr(t,:) = [0 1 0]; % green
                end
            else
                sz(t) = 1;
                cr(t,:) = [0.75 0.75 0.75]; % grey
            end
        end
    end
end

%W = table(tumortype,genex,geney,gi,gj,dotx,doty,sz,pval,odds);

%%

figure('position',[0 600 1200 500])
hold on

yyaxis left
hcom = scatter(dotx, doty, sz*75, cr, 'filled','MarkerEdgeColor','k');

xlim([-(ng+1),ng+1])
ylim([0,ng+1])

for i = 2:ng
    text(0,ng - i + 1,genelist{i},'horizontalalignment','center','fontsize',12)
end

for i = 1:(ng-1)
    ht = text(i + 0.4,ng - i + 0.4,genelist{i},'horizontalalignment','left','fontsize',12);
    ht.Rotation = 45;
end

for i = 1:(ng-1)
    ht = text(- i - 0.4,ng - i + 0.4,genelist{i},'horizontalalignment','right','fontsize',12);
    ht.Rotation = 360-45;
end

set(gca,'visible','off')
hold off
