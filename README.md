
# About SIREN in PET and CT

- Currently only the Reconstruction task in PET images doesn't work. 
- For all the other 3 cases the scripts work.

The main challenge with reconstruction in PET is that projection operations are required, but the output of the model must remain attached (i.e., connected to the computation graph) for gradients to flow properly during training. Therefore, we needed a library capable of performing projection while keeping the model output as an attached tensor.

We used Parallelproj for this purpose, as it supports differentiable projection operations. Note, however, that Parallelproj is only compatible with Unix-based systems — we used Oracle VM running Ubuntu 22.04 with 16 GB of RAM for our setup.

If the model output were detached, gradients would be lost and no learning would be possible.

We were not able to fully explore the library or understand all of its functionalities, but the authors provide several helpful resources. They specifically recommend reviewing the 2023 MIC short course they presented. One useful reference is the section on “forward and back projection” with custom linear operators (such as Parallelproj projectors), defined here:

🔗 https://github.com/gschramm/2023-MIC-ImageRecon-Shortcourse/blob/main/layers.py

Additionally, we included a PDF of our email exchange with the author, which may be helpful. It’s titled "Email to Georg_1.pdf" and can be found in the "Presentations, Docs and Results" folder.

## Documentation

[Documentation](https://linktodocumentation)

 - The folder "Presentations, Docs and Results" contains the experimentations we realized, with the quanlitative results in separeted, correctly named folders to the hyperparameters. The quanlitative results are in an excel.
 - Its a bit disorganized, but you can guide yourself by looking at the excel and the looking for the corresponding folder.
 - It also contations our some presentations and a first report.

 - THE MAIN FILES FOR TESTING ARE IN "Notebooks" folder :
    - CT denoising.ipynb for CT denoising
    - CT reconstruction.ipynb for CT reconstruction.
    - PET denoising.ipynb for PET denoising.
    - The other noteboooks in "Supplementary/Extra notebooks" are individual tests we did to to experiment more freely. You can look at them if you want, but they wont help a lot.

- For are first approach to the PET reconstruction task, the corresponding code is in the "CODE YOUNES" folder, into data, into phantome_after_augmentation, the notebook called "code_rec_TEP.ipynb".
    - The simulator_mmr_2d folder contains the algorithm to generate the data.

## Authors

- Younes Moussaoui
- Aureo Henrique E Silva Marques
- Mauricio Salim Gomez Chicre


