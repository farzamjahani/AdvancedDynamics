%%
function zdot=decoupled_equation(t,z)
q1=z(1);
q2=z(2);
q3=z(3);
q4=z(4);

qd1=z(5);
qd2=z(6);
qd3=z(7);
qd4=z(8);

A=[             15/4,  (3*sin(q3))/2,    (3*q2*cos(q3))/2,           cos(q4)/5;
    (3*sin(q3))/2,            3/2,                   0,      sin(q3 - q4)/5;
 (3*q2*cos(q3))/2,              0,          (3*q2^2)/2, (q2*cos(q3 - q4))/5;
        cos(q4)/5, sin(q3 - q4)/5, (q2*cos(q3 - q4))/5,                4/75];

B=[                                        0,                                        0,                                                   0,                                                                                       0;
                        (3*qd3*cos(q3))/2,                                        0, 3*q2*qd3 + (3*qd1*cos(q3))/2 + (qd4*cos(q3 - q4))/5,                                                                    (qd3*cos(q3 - q4))/5;
 (3*qd2*cos(q3))/2 - (3*q2*qd3*sin(q3))/2, (3*qd1*cos(q3))/2 + (qd4*cos(q3 - q4))/5,      -(q2*(15*qd1*sin(q3) + 2*qd4*sin(q3 - q4)))/10, (cos(q4)*(qd2*cos(q3) - q2*qd3*sin(q3)))/5 + (sin(q4)*(qd2*sin(q3) + q2*qd3*cos(q3)))/5;
                         -(qd4*sin(q4))/5,                    -(qd4*cos(q3 - q4))/5,                             (q2*qd4*sin(q3 - q4))/5,                        (q2*qd3*sin(q3 - q4))/5 - (qd2*cos(q3 - q4))/5 - (qd1*sin(q4))/5];
 dL_dq = [                                                                                         75 - 75*q1;
                               (2943*cos(q3))/200 - 50*q2 + (3*q2*qd3^2)/2 + (3*qd1*qd3*cos(q3))/2 + (qd3*qd4*cos(q3 - q4))/5 + 20;
 (3*qd1*qd2*cos(q3))/2 - (2943*q2*sin(q3))/200 + (qd2*qd4*cos(q3 - q4))/5 - (3*q2*qd1*qd3*sin(q3))/2 - (q2*qd3*qd4*sin(q3 - q4))/5;
                                  (q2*qd3*qd4*sin(q3 - q4))/5 - (qd1*qd4*sin(q4))/5 - (qd2*qd4*cos(q3 - q4))/5 - (981*sin(q4))/500];


qd=z(5:end);
qdd=A\(dL_dq - B*qd);

zdot=zeros(8,1);
zdot(1:4)=qd;
zdot(5:end)=qdd;



end