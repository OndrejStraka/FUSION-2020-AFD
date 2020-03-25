
PROCESS_DATA = 1;
if PROCESS_DATA
  disp('Processing data ...')
  load('decentralized-design-GPB1','nS','innerEdges','gamma','nodedec','V');
  V_GPB1{1} = reshape(gamma{1},nS{1});
  V_GPB1{2} = reshape(gamma{2},nS{2});
  load('decentralized-design-IMM','nS','innerEdges','gamma','nodedec','V');
  V_IMM{1} = reshape(gamma{1},nS{1});
  V_IMM{2} = reshape(gamma{2},nS{2});
  load('decentralized-design-GPB2','nS','innerEdges','gamma','nodedec','V');
  V_GPB2{1} = reshape(gamma{1},nS{1});
  V_GPB2{2} = reshape(gamma{2},nS{2});

  V_IMM_zip{1} = zeros(size(V_GPB1{1}));
  V_IMM_zip{2} = zeros(size(V_GPB1{2}));
  V_GPB2_zip{1} = zeros(size(V_GPB1{1}));
  V_GPB2_zip{2} = zeros(size(V_GPB1{2}));

  for im = 1:size(V_GPB1{1},1)/2
    for iv = 1:size(V_GPB1{1},2)
      for ip = 1:size(V_GPB1{1},3)
        V_IMM_zip{1}(im,iv,ip) = V_IMM{1}(im,im,iv,iv,ip);
        V_IMM_zip{2}(im,iv,ip) = V_IMM{2}(im,im,iv,iv,ip);
        V_GPB2_zip{1}(im,iv,ip) = V_GPB2{1}(im,im,iv,iv,ip);
        V_GPB2_zip{2}(im,iv,ip) = V_GPB2{2}(im,im,iv,iv,ip);
      end
    end
  end
  save test_Bellman_data
else
  disp('Loading data ...')
  load test_Bellman_data
end

figure
for ip = 1:size(V_GPB1{1},3)
  subplot(2,3,1)
  surf(V_GPB1{1}(:,:,ip))
  title({sprintf('%f',states{1}{3}(ip));'GPB1, subsystem 1'}) 
  subplot(2,3,2)
  surf(V_IMM_zip{1}(:,:,ip))
  title('IMM, subsystem 1') 
  subplot(2,3,3)
  surf(V_GPB2_zip{1}(:,:,ip))
  title('GPB2, subsystem 1') 
  subplot(2,3,4)
  surf(V_GPB1{2}(:,:,ip))
  title('GPB1, subsystem 2') 
  subplot(2,3,5)
  surf(V_IMM_zip{2}(:,:,ip))
  title('IMM, subsystem 2') 
  subplot(2,3,6)
  surf(V_GPB2_zip{2}(:,:,ip))
  title('GPB2, subsystem 2') 
  pause
end
