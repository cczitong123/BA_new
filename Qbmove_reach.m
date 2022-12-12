function [l, l_x, l_xx, l_u, l_uu, l_ux] = Qbmove_reach(x,u,t)
            %Rtf  = diag([w_tf, 0, 0, 0, 0, 0]);
          
            L = [0.45;0.45];
            x1 = x(1);
            x2 = x(2);
            xdot1 = x(3);
            xdot2 = x(4);
            l1 = L(1);
            l2 = L(2);
            l2_r = 0.45;
            g = 10;
            y0 = 1;
            p = 20;
            T = 2;


            x_m = l1*cos(x1) + l2*cos(x1+x2);
            xdot_m = -l1*sin(x1)*xdot1 - l2_r*sin(x1+x2)*(xdot1+xdot2);
            y_m = l1*sin(x1) + l2*sin(x1+x2);
            ydot_m = -l1*cos(x1)*xdot1 - l2_r*cos(x1+x2)*(xdot1+xdot2);
            T_m = (1/g)*((ydot_m) + sqrt(ydot_m^2 + 2*g*(y_m - y0)));
            
            d = x_m + xdot_m*T_m; %performance part
            
            %j = -d + 1/2*int((theta_dot).^2, t, T, 0);
            j = -d;

            %{
                Hf = [1,0,0,0;
                0,1,0,0;
                0,0,1,0;
                0,0,0,1];
            Hr = [1,1,0,0;
                1,1,0,0;
                0,0,1,1;
                0,0,1,1];
            %}
            if (isnan(u))
                % final cost
                %fh = @(x)obj.costf_reach(x);
                
                l = j; %+ Hf(5)*u(1)^2/2 + Hf(6)*u(2)^2/2 + Hf(7)*u(4)^2/2 + Hf(8)*u(5)^2/2;
                if nargout > 1
                    l_x = jacobian(j, x);
                    l_xx = hessian(j, x);
                    %l_x = (x')*Hf ;
                    %l_xx = Hf;
                end
            else
                % running cost
                %para = [];
                %para.w_t = w_t;
                %para.w_e = w_e;
                %fl = @(x,u,t) obj.costr_reach(x,u,t);
                %l = fl(x,u,t);
                l = j;
                if nargout > 1
                    %{
                    % finite difference
                    %flJ=@(x,u,t)J_cost_fd ( fl, x, u, t );
                    %[l_x ,l_u      ] = flJ ( x, u, t );
                    %flH =@(x,u,t)H_cost_fd  ( flJ, x, u, t );
                    %[l_xx,l_uu,l_ux] = flH  ( x, u, t );
                    
                    l_x = [0 0 0 0 0 0 0 0];
                    l_u = u'*Hr;
                    l_xx = zeros(4,4);
                    l_ux = zeros(4,4);
                    l_uu = Hr;
                    %}
                    l_x = get_jacobian_fd(j, x);
                    l_xx = get_hessian_fd(j, x);
                    l_u = get_jacobian_fd(j, u);
                    l_uu = get_hessian_fd(j, u);
                    l_ux = get_jacobian_fd(get_jacobian_fd(j, u), x);
                end
            end




end

