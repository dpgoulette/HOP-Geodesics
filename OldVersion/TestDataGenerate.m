x_1 = 3 + 1*randn(2000,1);
y_1 = 2 + 1*randn(2000,1);
x_2 = 4.5 + 5.5*rand(200,1);
y_2 = -.5.*x_2 + 4;
x_3 = 10 + .25*randn(200,1);
y_3 = -1 + .25*randn(200,1);

dust_x=-1+12*rand(100,1);
dust_y=-2+7*rand(100,1);

TestData = [x_1,y_1; x_2, y_2; x_3,y_3; dust_x,dust_y];

fprintf('\nWhat do you want to name the test data file (including the file extension)?\n')
namestring = input('Type the file name in single quotes: ');

save(namestring,'TestData','-ascii')
