robot = Qbmove2Dof();

x0 = zeros(4,1);
u0 = zeros(4,1);
%u0 = [pi/2;pi/2;pi/2;pi/2];
%u0 = [pi/3;pi/2;pi/2;pi/2];     %[q1_for motor1, q2_for motor1, q1_for motor2, q2_for motor2]
%robot.step(x0,ones(6,1),0.02)
u = repmat(u0,1,100);
%for i=1:100
%    u(1,i) = u(1,i)+i*pi/180;
%    u(2,i) = u(2,i)+2*i*pi/180;
%    u(3,i) = u(3,i)+i*pi/180;
%    u(4,i) = u(4,i)+2*i*pi/180;
%end

xsim = robot.simulate_feedforward(x0,u,0.02);

%%
%global mccpvd2_fh mccpvd2_ph1 mccpvd2_ph2
robot.plot(x0);

%%
robot.update(xsim(:,end));
%%

robot.plot([-pi/3; 2*pi/3])

%%
frames = robot.animate(xsim,0.02);
robot.savegif(frames, 'Qbmove_animation')
%%
movie(frames,1,50)

%%
gifname = 'Qbmove2_animation.gif';
for idx = 1:length(frames)
    im = frame2im(frames(idx));
    [A,map] = rgb2ind(im,256);
    if idx == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.02);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.02);
    end
end

