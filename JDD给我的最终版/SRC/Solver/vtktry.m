load mri  
 D = squeeze(D);  
 vtkwrite('mri', 'structured_points', 'mri', D)