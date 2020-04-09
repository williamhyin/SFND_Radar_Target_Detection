# Radar Target Generation and Detection

Project for Udacity's Sensor Fusion Engineer Nanodegree Program

### Project Goal

- Configure the FMCW waveform based on the system requirements.
- Define the range and velocity of target and simulate its displacement.
- For the same simulation loop process the transmit and receive signal to determine the beat signal
- Perform Range FFT on the received signal to determine the Range
- Towards the end, perform the CFAR processing on the output of 2nd FFT to display the target.

### Overview

![](https://williamhyin-1301408646.cos.ap-shanghai.myqcloud.com/img/20200409113409.png)

### Steps

1. Radar specifications

   ```matlab
   %% Radar Specifications 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Frequency of operation = 77GHz
   % Max Range = 200m
   % Range Resolution = 1 m
   % Max Velocity = 100 m/s
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   max_range=200;
   c = 3e8;
   range_resolution= 1;
   
   %Operating carrier frequency of Radar 
   fc= 77e9;             %carrier freq
   ```

2. Target specifications

   ```matlab
   %% User Defined Range and Velocity of target
   % *%TODO* :
   % define the target's initial position and velocity. Note : Velocity
   % remains contant
   target_pos=100;
   target_speed=30;
   ```

3. FMCW Waveform Generation

   In this project, we will designing a Radar based on the given system requirements (above).

   Max Range and Range Resolution will be considered here for waveform design.

   ```matlab
   % *%TODO* :
   %Design the FMCW waveform by giving the specs of each of its parameters.
   % Calculate the Bandwidth (B), Chirp Time (Tchirp) and Slope (slope) of the FMCW
   % chirp using the requirements above.
   
   B_sweep = c/(2*range_resolution); %Calculate the Bandwidth (B)
   
   T_chirp = 5.5*2*max_range/c;
   
   slope=B_sweep/T_chirp;
   ```

   Then, we need  simulate the signal propagation and move target scenario.

   ![](https://williamhyin-1301408646.cos.ap-shanghai.myqcloud.com/img/20200409114658.png)

   ```matlab
   %% Signal generation and Moving Target simulation
   % Running the radar scenario over the time. 
   
   for i=1:length(t)         
       
       
       % *%TODO* :
       %For each time stamp update the Range of the Target for constant velocity. 
       
       r_t(i) = target_pos+(target_speed*t(i));
       td(i) = 2*r_t(i)/c; % Time delay 
       
       % *%TODO* :
       %For each time sample we need update the transmitted and
       %received signal. 
       Tx(i) = cos(2 * pi * (fc * t(i) + slope * (t(i)^2)/2));
       Rx(i) = cos(2 * pi * (fc * (t(i) - td(i)) + slope * ((t(i)-td(i))^2)/2));
       
       % *%TODO* :
       %Now by mixing the Transmit and Receive generate the beat signal
       %This is done by element wise matrix multiplication of Transmit and
       %Receiver Signal
       Mix(i) = Tx(i) .* Rx(i);
       
   end
   ```

4. Range measurement

   The 1st FFT output for the target located at 100 meters

   ![](https://williamhyin-1301408646.cos.ap-shanghai.myqcloud.com/img/20200409115751.png)

   

5. Range and Doppler measurement

   2st FFT will generate a Range Doppler Map as seen in the image below and it will be given by variable ‘RDM’. 

   ![](https://williamhyin-1301408646.cos.ap-shanghai.myqcloud.com/img/20200409115935.png)

6.  CFAR implementation

   The 2D CFAR is similar to 1D CFAR, but is implemented in both dimensions of the range doppler block. The 2D CA-CFAR implementation involves the training cells occupying the cells surrounding the cell under test with a guard grid in between to prevent the impact of a target signal on the noise estimate.

   ![](https://williamhyin-1301408646.cos.ap-shanghai.myqcloud.com/img/20200409120132.png)

   1. Select the number of Training Cells and Guard Cells in both the dimensions and set offset of threshold

      ```matlab
      % *%TODO* :
      %Select the number of Training Cells in both the dimensions.
      
      Tr=10;
      Td=8;
      
      % *%TODO* :
      %Select the number of Guard Cells in both dimensions around the Cell under 
      %test (CUT) for accurate estimation
      
      Gr=4;
      Gd=4;
      
      % *%TODO* :
      % offset the threshold by SNR value in dB
      
      off_set=1.4;
      ```

   2. Slide Window through the complete Range Doppler Map

      ```matlab
      % *%TODO* :
      %design a loop such that it slides the CUT across range doppler map by
      %giving margins at the edges for Training and Guard Cells.
      %For every iteration sum the signal level within all the training
      %cells. To sum convert the value from logarithmic to linear using db2pow
      %function. Average the summed values for all of the training
      %cells used. After averaging convert it back to logarithimic using pow2db.
      %Further add the offset to it to determine the threshold. Next, compare the
      %signal under CUT with this threshold. If the CUT level > threshold assign
      %it a value of 1, else equate it to 0.
      
      
      % Use RDM[x,y] as the matrix from the output of 2D FFT for implementing
      % CFAR
      
      RDM = RDM/max(max(RDM)); % Normalizing
      
      % *%TODO* :
      % The process above will generate a thresholded block, which is smaller 
      %than the Range Doppler Map as the CUT cannot be located at the edges of
      %matrix. Hence,few cells will not be thresholded. To keep the map size same
      % set those values to 0. 
      
      %Slide the cell under test across the complete martix,to note: start point
      %Tr+Td+1 and Td+Gd+1
      for i = Tr+Gr+1:(Nr/2)-(Tr+Gr)
          for j = Td+Gd+1:(Nd)-(Td+Gd)
              %Create a vector to store noise_level for each iteration on training cells
              noise_level = zeros(1,1);
              %Step through each of bins and the surroundings of the CUT
              for p = i-(Tr+Gr) : i+(Tr+Gr)
                  for q = j-(Td+Gd) : j+(Td+Gd)
                      %Exclude the Guard cells and CUT cells
                      if (abs(i-p) > Gr || abs(j-q) > Gd)
                          %Convert db to power
                          noise_level = noise_level + db2pow(RDM(p,q));
                      end
                  end
              end
              
              %Calculate threshould from noise average then add the offset
              threshold = pow2db(noise_level/(2*(Td+Gd+1)*2*(Tr+Gr+1)-(Gr*Gd)-1));
              %Add the SNR to the threshold
              threshold = threshold + off_set;
              %Measure the signal in Cell Under Test(CUT) and compare against
              CUT = RDM(i,j);
              
              if (CUT < threshold)
                  RDM(i,j) = 0;
              else
                  RDM(i,j) = 1;
              end
              
          end
      end
      
      RDM(RDM~=0 & RDM~=1) = 0;
      ```

      The output of the 2D CFAR process，a peak and spread centered at 100m in range direction and 30 m/s in the doppler direction.

      ![](https://williamhyin-1301408646.cos.ap-shanghai.myqcloud.com/img/20200409121002.png)