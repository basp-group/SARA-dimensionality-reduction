function [u, v, Nm, uvidx, aW, nW] = util_gen_block_structure(uw, vw, aWw, nWw, param)

Nm = length(uw);

if param.use_manual_frequency_partitioning
    %% frequency clustering

    indexu = quantiz(uw, param.fpartition_y);
    indexv = quantiz(vw, param.fpartition_x);
    indexuv = length(param.fpartition_y)*indexu + indexv;
    indexuv = indexuv+1;

    pno = length(param.fpartition_x) * length(param.fpartition_y);

    u = cell(pno, 1);
    v = cell(pno, 1);
    
    uvidx = cell(pno, 1);
    uvidxw = 1:Nm;
    uvidxw = uvidxw(:);

    for q = 1:pno
        u{q} = uw(indexuv==q);
        v{q} = vw(indexuv==q);
        
        uvidx{q} = uvidxw(indexuv==q);
    end
end

if param.use_manual_partitioning
    %% data clustering

    pno = length(param.partition);

    u = cell(pno, 1);
    v = cell(pno, 1);
    
    Rp = 0;
    for q = 1:pno
        u{q} = uw(Rp+1:Rp+param.partition(q));
        v{q} = vw(Rp+1:Rp+param.partition(q));
        uvidx{q} = Rp+1:Rp+param.partition(q);
        uvidx{q} = uvidx{q}(:);
        Rp = Rp + param.partition(q);
    end
end

if param.use_uniform_partitioning
    pno = param.uniform_partitioning_no;
    
    Np = floor(Nm/pno) * ones(pno, 1);
    rand_loc = randi([1, pno], Nm - pno*floor(Nm/pno), 1);
    Np(rand_loc) = Np(rand_loc) + 1;
    
    u = cell(pno, 1);
    v = cell(pno, 1);
    
    uw_ = uw;
    vw_ = vw;
    
    uvidx = cell(pno, 1);
    uvidxw = 1:Nm;
    uvidxw = uvidxw(:);
    
    for q = 1:pno
        indexuv = [];
        while length(indexuv) < Np(q)
            indexuv = [indexuv; randi([1, length(uw_)], Np(q), 1)];
            indexuv = unique(indexuv);
        end
        indexuv = indexuv(1:Np(q));
        u{q} = uw_(indexuv);
        v{q} = vw_(indexuv);
        uvidx{q} = uvidxw(indexuv);
        uw_(indexuv) = [];
        vw_(indexuv) = [];
        uvidxw(indexuv) = [];
    end
end

if param.use_equal_partitioning
    pno = param.equal_partitioning_no;
    fpno = util_divisor(pno);
    xno = fpno(find(fpno >= sqrt(pno), 1));
    yno = pno/xno;
    

    su = ceil(yno*tiedrank(uw)/length(uw));
    
    
    
    u = cell(pno, 1);
    v = cell(pno, 1);

    uvidx = cell(pno, 1);
    uvidxw = 1:Nm;
    uvidxw = uvidxw(:);
    
    for q = 1:yno
        u_ = uw(su == q);
        v_ = vw(su == q);
        uvidx_ = uvidxw(su == q);
        
        sv = ceil(xno*tiedrank(v_)/length(v_));
        for k = 1:xno
            u{xno * (q - 1) + k} = u_(sv == k);
            v{xno * (q - 1) + k} = v_(sv == k);
            uvidx{xno * (q - 1) + k} = uvidx_(sv == k);
        end
    end
end

if param.use_density_partitioning
    pno = param.density_partitioning_no;
    
    Np = floor(Nm/pno) * ones(pno, 1);
    rand_loc = randi([1, pno], Nm - pno*floor(Nm/pno), 1);
    Np(rand_loc) = Np(rand_loc) + 1;
    
    [~, sidx] = sort(aWw);
    
    u = cell(pno, 1);
    v = cell(pno, 1);
    
    Rp = 0;
    for q = 1:pno
        u{q} = uw(sidx(Rp+1:Rp+Np(q)));
        v{q} = vw(sidx(Rp+1:Rp+Np(q)));
        uvidx{q} = sidx(Rp+1:Rp+Np(q));
        uvidx{q} = uvidx{q}(:);
        Rp = Rp + Np(q);
    end
end

aW = cell(pno, 1);
for k = 1:pno
    aW{k} = aWw(uvidx{k});
end

nW = cell(pno, 1);
for k = 1:pno
    nW{k} = nWw(uvidx{k});
end

end