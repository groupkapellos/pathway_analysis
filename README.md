## **Pathway Gene Finder Shiny App**

This is a Shiny application for finding pathways based on selected genes and viewing associated gene counts through a bar plot. The app allows users to explore pathways related to specific genes and visualize the number of genes associated with each pathway.

### **Prerequisites**

To run this Shiny app on your local machine, you will need the following software and packages installed:

**R:** Make sure you have R installed. If not, download and install it from CRAN.

**RStudio:** Although not strictly necessary, RStudio is a great integrated development environment (IDE) for R. You can download it from RStudio.

**Git:** You will also need to have Git installed for cloning the repository. If you're new to Git or prefer a simpler method, we recommend using GitHub Desktop.

### **Installing Git** 

**Using git:** Download and install Git from [here](https://git-scm.com/). Once installed, follow the installing dependencies section below to clone the repository.

**Using GitHub Desktop:** GitHub Desktop is a simpler alternative for users who are not familiar with git.

1.  Download and install GitHub Desktop from [here](https://desktop.github.com/).
2.  Once installed, open GitHub Desktop and sign in with your GitHub account.
3.  Go to File \> Clone repository.
4.  In the URL tab, paste the repository URL: `https://github.com/groupkapellos/pathway_analysis.git`
5.  Choose the folder where you want to clone the repository (e.g., Documents) and click Clone.

### **Installing Dependencies** 

This app requires several R packages. To make the installation process easier, we've provided a script that will automatically install the required packages for you.

**Step 1:** Download the Repository First, clone or download this repository to the folder where you want to store the files.

If you have Git installed, you can clone the repository by running:

```{bash}
git clone https://github.com/groupkapellos/pathway_analysis.git
```

**Step 2:** Install Required Packages Navigate to the folder where you downloaded the repository (for example, pathway-gene-finder/). Inside the repository, you'll find the *install_dependencies.R* script.

1.  Run the script to install the required R packages:

2.  Open R or RStudio.

3.  Set the working directory to the folder containing the repository. You can do this with the following command in R:

4.  run the command:

    ```{r}
    setwd("path/to/pathway_project") 
    ```

5.  Run the install_dependencies.R script:

    ```{r}
    source("install_dependencies.R")
    ```

This will check if the required packages are already installed. If any packages are missing, the script will install them for you.

**Step 3:** Launch the App Once the dependencies are installed, you can run the Shiny app:

In the same R session or RStudio, run the following command:

```{r}
shiny::runApp("path/to/shiny_app.R")
```

Alternatively, you can open and run the script run_pathway_app.R. Either of these methods will launch the app in your default web browser.
