clc; clear; close all;

% Data Generation
n = [10 50 100 300 500];
m = [5 30 80 250 450];

EndCond = 10^-8;
MaxIter = 10^6;
     
k = zeros(size(n,2),size(m,2));
f1=figure;
for i=1:size(n,2)
    for j=1:size(m,2)
        A = rand(m(j),n(i));
        x_seed = rand(n(i),1);% Guarantee Feasibility
        b = A*x_seed;
        
        % b = rand(m,1);

        % cvx_begin quiet
        %     variables x(n)
        %     minimize 0
        %     subject to
        %         A*x==b;
        %         x>=0;
        % cvx_end

        % if strcmp('Infeasible',cvx_status) == 1
        %     fprintf("Problem Infeasible\n");
        %     return
        % end
        
        P_S1 = @(x) x - A'*((A*A')\(A*x-b));
        P_S2 = @(x) max(x,0);

        %Alternating projection algorithm
        
        %Initialization
        x_0 = rand(n(i),1)*50;
        x_k = x_0;
        
        while(1)
            if (k(i,j) > MaxIter)
                break;
            end

            k(i,j) = k(i,j)+1;
            x_k = P_S2(P_S1(x_k));
            if max(A*x_k-b)<EndCond && min(x_k)>0
                break;
            end
        end
    end
    semilogy(m,k(i,:),'-*','DisplayName',strcat('n=',num2str(n(i))) );
    hold on;
    pause(0.01);
end
xlabel('\fontsize{18} m','interpreter','tex');
ylabel('\fontsize{18} Iterations','interpreter','tex');
legend show

saveas(f1,fullfile('D:\Documents\Tuc\HMMY\10th Semester\ConvexOptimization','P2-fig2.png'));