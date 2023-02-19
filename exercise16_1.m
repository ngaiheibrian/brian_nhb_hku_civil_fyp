clear all
close all

%Parameter kN/cm/s
mass=0.2591;
mass_m=mass*eye(5);
stiff=100; 
stiff_m=stiff*[ 2 -1 0 0 0;
               -1 2 -1 0 0;
               0 -1 2 -1 0;
               0 0 -1 2 -1;
               0 0 0 -1 1]

modal_damp_ratio=0.05;
timestep=0.01;
duration=20;

%Linear_acceleration
newmark_y=1/2; 
newmark_beta=1/6;

total_step= duration/timestep;



%Initital Calculations_parameters
%--Eigenvalue problem-->find natural frequencies and mode
[fq_sq, mode]=find_eign(mass_m,stiff_m);
fq=sqrt(fq_sq);
disp("Natural frequency:");
disp(fq);
disp("Mode:");
disp(mode);



%--Calculate the modal mass, modal stiffness, modal damp
modal_mass= transpose(mode)*mass_m*mode;

modal_stiffness= transpose(mode)*stiff_m*mode;


modal_damp=modal_damp_ratio*2*modal_mass*fq;

disp("modal_mass")
disp(modal_mass);
disp("modal_stiffness");
disp(modal_stiffness);
disp("modal_damp");
disp(modal_damp);



disp("force coe:")
force_coe=transpose(mode)*(-mass_m)*[1; 1; 1; 1; 1];
disp(force_coe);


%Initital Calculations
a1=(1/(newmark_beta*(timestep^2)))*modal_mass+(newmark_y/(newmark_beta*timestep))*modal_damp;
a2=(1/(newmark_beta*(timestep)))*modal_mass+((newmark_y/newmark_beta)-1)*modal_damp;
a3=((1/(2*(newmark_beta)))-1)*modal_mass+timestep*((newmark_y/(2*newmark_beta))-1)*modal_damp;

K_h=modal_stiffness+a1;



time=0;
p=zeros(5,1);

u=zeros(5,1);

q=zeros(5,1);
q_v=zeros(5,1);
q_a=zeros(5,1);

q_new=zeros(5,1);
q_v_new=zeros(5,1);
q_a_new=zeros(5,1);




%Calculations for each time step i

for i= 0:total_step
    
time= time+ timestep;
disp("time:");
disp(time);

p=force_coe*find_acce(time)+a1*q+a2*q_v+a3*q_a;
q_new=inv(K_h)*p;



q_v_new=(newmark_y/(newmark_beta*timestep))*(q_new-q)+(1-(newmark_y/newmark_beta))*q_v+timestep*(1-newmark_y/(2*newmark_beta))*q_a;
q_a_new=(1/(newmark_beta*(timestep^2)))*(q_new-q)-(1/(newmark_beta*timestep))*q_v-(1/(2*newmark_beta)-1)*q_a;


q=q_new;
q_v=q_v_new;
q_a=q_a_new;

u=mode*q;


disp(q);





%store data

data_q(i+1,1:5)=q;
data_u(i+1,1:5)=u;
data_time(i+1)=time;
data_force(i+1, 1:5)=force_coe*find_acce(time);










end

%Ploting

hold on
plot(data_time,data_u);
xlabel('Time per second');
ylabel('Displacement(cm)');
title('The displacement of the structure over time');
legend('X1','X2','X3','X4','X5');
hold off

figure;

hold on
plot(data_time,data_q);
xlabel('Time per second');
ylabel('Q');
title('The Q of the structure over time');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5');
hold off

figure;
hold on
plot(data_time,data_force);
xlabel('Time per second');
ylabel('Force');
title('The Force of the structure over time');
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5');
hold off

%Table making & Export to Excel
T= table(data_time(:),data_q(:,1),data_q(:,2),data_u(:,1),data_u(:,2),data_u(:,3),data_u(:,4),data_u(:,5));

T.Properties.VariableNames([1 2 3 4 5 6 7 8])= {'Time','q1','q2','x1','x2','x3','x4','x5'};

filename='behaviour2.xlsx';

writetable(T,filename,'Sheet',1,'Range','A1');
