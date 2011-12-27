load A;
B = A>25; % Randomly Chosen Threshold
A = A./max(max(max((A))));

for i = 1:size(B,3)
    imagesc([A(:,:,i) B(:,:,i)]);
    axis image;
    title('Orignial PET Image vs Threshold');
    drawnow;
end;

Extract = FloodFill3D(B,60);

for i = 1:size(B,3)
    imagesc([Extract(:,:,i) B(:,:,i)]);
    title('Extracted Image vs Threshold');
    drawnow;
end;