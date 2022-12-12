classdef Qbmove_cost
    
    properties

    end
    methods
    end

    methods(Static)


        function [l_x] = J_cost_final(l,x)

            dimX = size(x,1);
            delta=1e-6;
            
            l_x  = zeros(dimX,1);
            for i=1:dimX
	            dx=zeros(dimX,1); 
                dx(i)=delta;
            
	            lxp = l( x+dx );
	            lxm = l( x-dx );
            
	            l_x(i) = (lxp-lxm)/(2*delta);
            end

        end

        function [l_xx] = H_cost_final(J, x)


            dimX = size(x,1);
            
            delta=1e-6;
            
            l_xx = zeros(dimX,dimX);

            for i=1:dimX
	            dx=zeros(dimX,1); 
                dx(i)=delta;
            
	            lxxp = J( x+dx );
	            lxxm = J( x-dx );
            
	            l_xx(:,i) = (lxxp-lxxm)/(2*delta);
            end

        end
 

         function [l_x,l_u] = J_cost_fd ( l, x, u, t )

            dimX = size(x,1);
            dimU = size(u,1);
            
            delta=1e-6;
            
            l_x  = zeros(dimX,1);
            for i=1:dimX
	            dx=zeros(dimX,1); dx(i)=delta;
            
	            lxp = l( x+dx, u, t );
	            lxm = l( x-dx, u, t );
            
	            l_x(i) = (lxp-lxm)/(2*delta);
            end
%             dx=ones(dimX,1)*delta; 
% 
%             
% 	        lxp = l( x+dx, u, t );
% 	        lxm = l( x-dx, u, t );
%             
% 	        l_x = (lxp-lxm)/(2*delta);
%             l_u  = zeros(dimU,1);
            for i=1:dimU
	            du=zeros(dimU,1); du(i)=delta;
            
	            lup = l( x, u+du, t );
	            lum = l( x, u-du, t );
            
	            l_u(i) = (lup-lum)/(2*delta);
            end
% 	         du=ones(dimU,1)*delta;
%             
% 	         lup = l( x, u+du, t );
% 	         lum = l( x, u-du, t );
%             
% 	         l_u = (lup-lum)/(2*delta);      

         end

         function [l_xx,l_uu,l_ux] = H_cost_fd (J, x, u, t )

            dimX = size(x,1);
            dimU = size(u,1);
            
            delta=1e-6;
            
            l_xx = zeros(dimX,dimX);
            l_ux = zeros(dimU,dimX);
            for i=1:dimX
	            dx=zeros(dimX,1); dx(i)=delta;
            
	            [lxxp,luxp] = J( x+dx, u, t );
	            [lxxm,luxm] = J( x-dx, u, t );
            
	            l_xx(:,i) = (lxxp-lxxm)/(2*delta);
	            l_ux(:,i) = (luxp-luxm)/(2*delta);
            end
%             dx=ones(dimX,1)*delta; 
% 
%            
% 	        [lxxp,luxp] = J( x+dx, u, t );
% 	        [lxxm,luxm] = J( x-dx, u, t );
%             
% 	        l_xx = (lxxp-lxxm)/(2*delta);
% 	        l_ux = (luxp-luxm)/(2*delta);
%             l_uu = zeros(dimU,dimU);
            for i=1:dimU
	            du=zeros(dimU,1); du(i)=delta;
            
	            [dummy,luup] = J( x, u+du, t );
	            [dummy,luum] = J( x, u-du, t );
            
	            l_uu(:,i) = (luup-luum)/(2*delta);
            end

%             du=ones(dimU,1)*delta;
%             [dummy,luup] = J( x, u+du, t );
% 	        [dummy,luum] = J( x, u-du, t );
%             
% 	        l_uu = (luup-luum)/(2*delta); 
         end

         
         function c = costf_reach(x)
             
            dt = 0.2;
            w_t = 1e3;
            w_tf = w_t*dt;
            d = distance(x);
            c = -d*w_tf ;
            
         end


         function j = costr_reach(x,u,t)

            w_t = 1e3;
            w_e = 0;
%             w_e = 1e0;

            %energy part
%             ep = [(u(1) + u(2))/2;(u(1) + u(2))/2;(u(3) + u(4))/2;(u(3) + u(4))/2];
            ep = [(u(1) + u(2))/2;0;(u(3) + u(4))/2;0];
            x_state  = [x(1);0;x(2);0];
%             j_energy = (u - ep)'*(u - ep);
            j_energy = (x_state - ep)'*(x_state - ep);
            %j = -d + 1/2*int((u-ep)'*(u-ep), t, 2, 0);
            j = w_e*j_energy;

         end



         function [l, l_x, l_xx, l_u, l_uu, l_ux] = Qbmove_reach(x,u,t)
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);

            if (isnan(u))
                % final cost
                fh = @(x)Qbmove_cost.costf_reach(x);
                l = fh(x);
                
            
                if nargout > 1


                    J_cost = @( x)Qbmove_cost.J_cost_final(fh, x);
                    l_x = J_cost(x);

                    H_cost = @( x, u, t)Qbmove_cost.H_cost_final(J_cost, x);
                    
                    l_xx = H_cost(x);

                    %l_x = (x')*Hf ;
                    %l_xx = Hf;
                end
            else
                % running cost
                %para = [];
                %para.w_t = w_t;
                %para.w_e = w_e;
                fl = @(x,u,t) Qbmove_cost.costr_reach(x,u,t);
                l = fl(x,u,t);
       
                if nargout > 1
                   
                    J_cost = @( x, u, t)Qbmove_cost.J_cost_fd(fl, x, u, t);
                    [l_x, l_u] = J_cost( x, u, t);

                    H_cost = @( x, u, t)Qbmove_cost.H_cost_fd(J_cost, x, u, t);
                    
                    [l_xx,l_uu,l_ux] = H_cost(x, u, t);
%                     l_x = J_cost(1);
%                     l_u = J_cost(2);
%                     l_xx = H_cost(1);
%                     l_uu = H_cost(2);
%                     l_ux = H_cost(3);

                    
%                     l_x = [0 0 0 0 0 0 0 0];
%                     l_u = u'*Hr;
%                     l_xx = zeros(4,4);
%                     l_ux = zeros(4,4);
%                     l_uu = Hr;

                   
                   
                end
            end

         end
    end 
end

function d = distance(x)

            x1 = x(1);
            x2 = x(2);
            xdot1 = x(3);
            xdot2 = x(4);

            L = [0.3;0.3];
            l1 = L(1);
            l2 = L(2);
            g = 10;

            y0 = 1;

            x_m = l1*cos(x1) + l2*cos(x1+x2);
            xdot_m = -l1*sin(x1)*xdot1 - l2*sin(x1+x2)*(xdot1+xdot2);
            y_m = l1*sin(x1) + l2*sin(x1+x2);
            ydot_m = l1*cos(x1)*xdot1 + l2*cos(x1+x2)*(xdot1+xdot2);
%             T_m = (1/g)*((ydot_m) + sqrt(ydot_m^2 + 2*g*(y_m + y0)));
            T_m = (1/g)*((ydot_m) + sqrt(ydot_m^2 + 2*g*(y_m - y0)));
            
            d = x_m + xdot_m*T_m; %performance part   
end

