%{
hw4_p2.m plots first and second order transfer functions.

@Author: Carlos Carrasquillo
@Date: June 20, 2019 09:30:00PM 
%}

clear; clc;                      % clears workspace and command window
syms s;                          % creates a symbolic variable

%% Initializes each Transfer Function in a Structure
% creates a new MATLAB structure to hold transfer functions
problem = 1;

switch problem  
  case 1                             % transfer functions for problem 2
    trans_func.a = tf(2, [1 1]);	 % defines the first transfer function
    trans_func.b = tf(-5, [1 5]);    % defines the second transfer function
    trans_func.c = tf(6, [6 2]);     % defines the third transfer function
  case 3                             % transfer functions for problem 1
    trans_func.a = tf(2, [1 2 1]);	 % defines the first transfer function
    trans_func.b = tf(-7, [1 7 7]);  % defines the second transfer function
    trans_func.c = tf(1, [3 6 9]);   % defines the third transfer function  
  otherwise 
    warning('Only problems 1 and 3 are hard-coded.');  % error statement
    return         % returns program counter to location of function call
end

%% Creates a New Figure 
fig = figure('name', 'Homework 4 Problem 2', ...% creates new figure object
       'position', [0 0 900 600]);     % with a specified size (900 x 600)
movegui('center')                      % places figure in center of display
sgtitle('Step Responses of Transfer Functions'); % assigns a title to plots
cursor_mode = datacursormode(fig);     % returns data cursor mode object
set(cursor_mode, 'updatefcn', @set_data_cursor); % creates custom data tips
fig_pos = get(fig, 'position');        % fetches current position of figure

%% Plots Each Transfer Function
fields = fieldnames(trans_func); % cell aray of trans_func fields
n = length(fields);
tf_subplot = nan(1, length(fields));

for ii = 1:n                     % iterates for every trans func
  
  % determines the order of the transfer function 
  order = length(trans_func.(fields{ii}).Numerator{1}) - 1;

  if order ~= 1 && order ~= 2
    warning('This program is not optimized for orders higher than 2.');
    return           % returns program counter to location of function call
  end
  
  % lines 16-19 create the LaTeX title for each transfer function
  num = poly2sym(trans_func.(fields{ii}).Numerator, s);   % format numer.
  den = poly2sym(trans_func.(fields{ii}).Denominator, s); % format denom.
  tf_sym = num / den;                                     % format num/den
  tf_title = latex(tf_sym);               % creates LaTex title
    
  tf_subplot(ii) = subplot(n, 1, ii);     % subplot for each tf
  [Y, T] = step(trans_func.(fields{ii})); % computes step response 
  tf_plot = plot(T, Y);                         % plots step response
  
  switch order   % selects based on the order of the system
    case 1       % displays time constant if the tf is a first-order system
      Kp = evalfr(trans_func.(fields{ii}), 0);  % calculates gain via G(0)
      data_tip = cursor_mode.createDatatip(tf_plot);    % creates data tip 
      
      sixty_three = (1 - 1/exp(1))*Kp;        % actual time constant value
      [~,index] = min(abs(Y - sixty_three));  % finds corresponding index
      tip_pos_Y = Y(index);                   % defines data tip Y-position
      tip_pos_T = T(index);                   % defines data tip T-position
      
      set(data_tip, 'position', [tip_pos_T tip_pos_Y]); % sets data tip
      
    % displays natural frequency (omega), damping ratio (zeta), and DC gain
    % (Kp) if the tf is a second-order system 
    case 2   
      %arrays containing the numerator and the denominator of the tf
      numer_coeff = trans_func.(fields{ii}).Numerator{1};   
      denom_coeff = trans_func.(fields{ii}).Denominator{1}; 
      
      %arrays of the numers and denoms scaled by the highest order term
      scaled_numer_coeff = numer_coeff / denom_coeff(1);
      scaled_denom_coeff = denom_coeff / denom_coeff(1); 
      
      omega = sqrt(scaled_denom_coeff(end));       % natural frequency calc
      Kp = scaled_numer_coeff(end) / omega^2;      % DC gain calc
      zeta = scaled_denom_coeff(2) / (2 * omega);  % damping ratio calc
  end
  
  subplot_pos = get(tf_subplot(ii), 'position');    % gets subplot position
  note_loc_x = subplot_pos(1) + subplot_pos(3) - .1;% textbox x-location
  note_loc_y = subplot_pos(2);                      % textbox y-location
  txt = annotation('textbox', [note_loc_x note_loc_y .1 .1], ...
      'FitBoxToText','on');                         % new textbox object
  txt.FontSize = 8;                                 % set font-size
  txt.EdgeColor = 'none';                           % no edge-color
  
  switch order  % selects based on the order of the system
    case 1
      set(txt, 'String', ['K_p = ' num2str(Kp)]);      % sets text in box
  
    case 2
      set(txt, 'String', ['K_p = ' num2str(Kp)]); 
  end
  
  %assigns the corresponding title to each plot
  title(sprintf('G(s) = $$ %s $$', tf_title), 'interpreter','latex');
end

%% Assigns Axes Labels
pos_1 = get(tf_subplot(1), 'position');   % gets position of first subplot
pos_n = get(tf_subplot(n), 'position');   % gets position of last subplot

height = pos_1(2) - pos_n(2) + (pos_1(4) * 2 / n); % computes new axes pos.

% creates a new axes object
axes('position',[pos_n(1) pos_n(2) pos_n(3) height],'visible','off');

ylabel('Amplitude','visible','on');       % sets y-axis label
xlabel(tf_subplot(n), 'Time');            % sets x-axis label

%% Sets Figure Caption for First-Order Transfer Functions
% creates a new axes object
axes('position',[pos_n(1) (pos_n(2)-.045) pos_n(3) 1],'visible','off');
caption = ['Note: The data-tips indicate the time constant' ...
          ' of a BIBO-stable system (63% of K_p).'];    % defines caption
set(gca,'FontSize', 8);                     % sets axis font-size
xlabel(caption, 'visible','on');            % sets x-axis label

%% Edit Cursor Axes Names
function labels = set_data_cursor(~,cursor)
%{
labels.m displays the data cursor position in a data tip.
 
@Param: 
    - cursor: the handle for the cursor_mode object

@Return: 
	- label: a cell array of character vectors containing the text to be 
             displayed with the tip text, 

@Author: Carlos Carrasquillo
@Date: June 20, 2019 11:00:00PM 
%}

pos = cursor.Position;      % determines the position of the cursor
% defines the 'x' and 'y' names for the data tip
labels = {['t',format(pos(1))], ['A',format(pos(2))]};
end

function formatted = format(to_format)
%{
format.m formats the data tip 'x' and 'y' values using TeX.

@Param:
    - to_format: the number to be formated

@Return: 
    - formatted: the formatted number with four decimal places

@Author: Carlos Carrasquillo
@Date: June 20, 2019 11:10:00PM 
%}
formatted = [' \color[rgb]{0 0.5 1}\bf' num2str(to_format,4) ...
    '\color[rgb]{.15 .15 .15}\rm']; %formats the number to four decimals
end

