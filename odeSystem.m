function ddy = odeSystem(~,y,reactionParameters)

ddy=zeros(3,1);

B = y(1); 
D = y(2);
C1 = y(3);
C2 = y(4);
%ModelParameters = [alpha, I, gamma, delta, epsilon, phi ];

v_mu_D = reactionParameters(1);
CC_D = reactionParameters(2);
v_delta_D = reactionParameters(3);
v_mu_B = reactionParameters(4);
CC_B = reactionParameters(5);
K_doc1 = reactionParameters(6);
v_delta_B = reactionParameters(7);
lambda = reactionParameters(8);
K_doc2 = reactionParameters(9);
delta_doc1 = reactionParameters(10);
delta_doc2 = reactionParameters(11);
rAB = reactionParameters(12);

mu_D = v_mu_D * (1-D/CC_D);
delta_D = v_delta_D * (1/(1+mu_D)); 
mu_B = v_mu_B * (1-B/CC_B)*C1/(C1+K_doc1) + v_mu_B * (1-B/CC_B)*C2/(C2+K_doc2);
delta_B = v_delta_B *(1/(1+mu_B));

lambda_v = lambda ;
 
ddy(1)= -delta_B*B + mu_B*B*C1 + mu_B*B*C2 - rAB*B*D/CC_D;
ddy(2)= -lambda_v*D - delta_D*D + mu_D*D;
ddy(3)= lambda_v*D  - mu_B*B*C1 - delta_doc1*C1;
ddy(4)= delta_D*D - mu_B*B*C2 - delta_doc2*C2;

