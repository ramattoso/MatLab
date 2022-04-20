function [XData, YData] = xydata(nnodes, A,B,C,D)
    xsize = B-A;
    ysize = D-C;
    if xsize == ysize || xsize>ysize
        N = nnodes; K = 1:nnodes; MMC = K(rem(N,K)==0); auxMMC = MMC;
        if size(MMC,2) < 3
            N = N + 1;
            K = 1:N; MMC = K(rem(N,K)==0);
        end
        m = MMC(floor((size(MMC,2)/2))+1);
        n = MMC(ceil((size(MMC,2)/2)));
        x = A:xsize/(m):B;
        y = C:ysize/n:D;
        aux=0;
        for i=1:size(x,2)-1
            for j=1:size(y,2)-1
                XData(j+aux,1) = A + xsize/(2*m)+(i-1)*xsize/m;
                YData(j+aux,1) = C + ysize/(2*n)+(j-1)*ysize/n;
            end
            aux=aux+j;
        end
        if size(auxMMC,2) < 3
            XData(aux-1) = [];
            YData(aux-1) = [];
        end
    else
        N = nnodes; K = 1:nnodes; MMC = K(rem(N,K)==0); auxMMC = MMC;
        if size(MMC,2) < 3
            N = N + 1;
            K = 1:N; MMC = K(rem(N,K)==0);
        end
        n = MMC(floor((size(MMC,2)/2))+1);
        m = MMC(ceil((size(MMC,2)/2)));
        x = A:xsize/m:B;
        y = C:ysize/(n):D;
        aux=0;
        for i=1:size(x,2)-1
            for j=1:size(y,2)-1
                XData(j+aux,1) = A + xsize/(2*m)+(i-1)*xsize/m;
                YData(j+aux,1) = C + ysize/(2*n)+(j-1)*ysize/n;
            end
            aux=aux+j;
        end
        if size(auxMMC,2) < 3
            XData(aux-1) = [];
            YData(aux-1) = [];
        end
    end
end