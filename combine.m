function Yp = combine(P,BPAs)
Yp = zeros(size(BPAs,1),1);
for i=1:size(BPAs,1)
    final_mass = [0 P.maxP 1];
    for j=1:size(BPAs,2)
        help_mass = zeros(size(final_mass,1)*size(BPAs{i,j}.sets,1),3);
        for i1=1:size(final_mass,1)
            for i2=1:size(BPAs{i,j}.sets,1)
                lb = max(final_mass(i1,1),BPAs{i,j}.sets(i2,1));
                ub = min(final_mass(i1,2),BPAs{i,j}.sets(i2,2));
                if lb >= ub
                    lb = 0;
                    ub = P.maxP;
                end
                help_mass(i2+(i1-1)*size(BPAs{i,j}.sets,1),:) = [lb ub final_mass(i1,3)*BPAs{i,j}.masses(i2)];
            end
        end
        r = unique(help_mass(:,1:2),'rows');
        final_mass = r;
        for i1=1:size(r,1)
            final_mass(i1,3) = sum(help_mass(ismember(help_mass(:,1:2),r(i1,:),'rows'),3));
        end
    end
    x = linspace(0,P.maxP,100);
    y = final_mass(:,3)'*(final_mass(:,1) < x & final_mass(:,2) > x);
    Yp(i) = x(find(cumsum(y == max(y))>=sum(y == max(y))/2,1));
end