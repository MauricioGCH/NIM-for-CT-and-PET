
# About SIREN in PET and CT

- Currently only the Reconstruction task in PET images doesn't work. 
- For all the other 3 cases the scripts work.

The problem with Reconstruction in PET is that the projection is need will the output of the model is still attach, an as such, a library that did that projection while the output was an attach tensor to the model was needed. parallelproj has that capability (Only compatibel with Unix systems, we used the Oracle VM with Ubuntu 22.04 setting 16gb of RAM)

If detached from the model, it would lose its gradient and no learning was possible. We weren't able to explore the totality of the library nor understand alls of its functionalities, but the authors have several resources to help undestarnd it.

They recommend to have a look at the 2023 MIC shortcourse that they gave. Specifically at the “forward and back projection” with custom linear operator (e.g. parallelproj projectors) as defined here:

 https://github.com/gschramm/2023-MIC-ImageRecon-Shortcourse/blob/main/layers.py

(I also added a pdf of our exchange, it might be helpful. It's called "Email to Georg_1.pdf" in the "Presentations, Docs and Results folder")

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


