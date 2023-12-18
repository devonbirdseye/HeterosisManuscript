#author: Devon Birdseye
#date: 11/16/2020

library(shiny)
library(ggplot2)
tmt.SL <- read.csv("tmt.SL.csv", stringsAsFactors = F)
cpm.SL <- read.csv("cpm.SL.csv", stringsAsFactors = F)
tmt.LB <- read.csv("tmt.LB.csv", stringsAsFactors = F)
cpm.LB <- read.csv("cpm.LB.csv", stringsAsFactors = F)
tmt.RIL <- read.csv("tmt.RIL.csv", stringsAsFactors = F)
cpm.RIL <- read.csv("cpm.RIL.csv", stringsAsFactors = F)
tmt.6H <- read.csv("tmt.6H.csv", stringsAsFactors = F)
cpm.6H <- read.csv("cpm.6H.csv", stringsAsFactors = F)

ScatterplotFun <- function(df, xcol, ycol, xlabel, ylabel, lims, Gene, title){
    #set limits
    df[[xcol]][df[[xcol]] > lims] <- lims
    df[[xcol]][df[[xcol]] < -lims] <- -lims
    df[[ycol]][df[[ycol]] > lims] <- lims
    df[[ycol]][df[[ycol]] < -lims] <- -lims
    #Set GOIs
    df$Category <- "Other"
    df[(df$Gene %in% Gene),"Category"] <- "Gene"
    #factor
    df$Category <- factor(df$Category, levels = c("Other", "Gene"))
    df <- df[(order(df$Category)),]
    #scatterplot
    df.scat <- ggplot(df, mapping=aes(x=get(xcol), y=get(ycol), color=Category, shape=Category))+
        geom_point(aes(color=Category), size=3.5)+
        scale_color_manual(values=c("black", "red"))+
        scale_shape_manual(values=c(1,19))+
        theme_classic()+
        xlim(c(-lims,lims))+
        ylim(c(-lims,lims))+
        geom_vline(xintercept=0, linetype="dashed", color = "gray60", size=0.5)+
        geom_hline(yintercept=0, linetype="dashed", color = "gray60", size=0.5)+
        xlab(paste0("log2 ", xlabel))+
        ylab(paste0("log2 ", ylabel))+
        theme(aspect.ratio=1)+
        coord_fixed(expand = F)+
        ggtitle(title)+
        theme(axis.text=element_text(size=10),
              legend.position="none", plot.title = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title =  element_text(size = 12))
    df.scat
}
#Define UI for application
ui <- fluidPage(
    titlePanel("Expression Visualization", windowTitle="Expression Visualization"),
    h4("Data published in Birdseye et al., 2021"),
    h5("Visualize expression heterosis for your gene(s) of interest"),
    h6("(If no points are plotted, the gene was not detected in the dataset)"),
    #Sidebar with a text input for gene and a button choice for "RNA" or "Protein"
    fluidRow(
        column(width = 6,
            wellPanel(
                textInput(inputId = "Gene", label = "Zea mays v4 gene accessions (e.g. Zm00001d002325) separated by space or comma", placeholder = "Zm00001d002325"),
                radioButtons(inputId = "dataset", label = "Dataset",
                         choices = list("Protein" = 1, "RNA" = 2), 
                         selected = 1)
            ),
        h5("Pearson correlation between expression heterosis and plant height heterosis from 6 hybrids and 8 RIL hybrids"),
        h6("(If no data appears, the gene was not detected in these datasets)"),
        #Correlation table
        tableOutput("Cor"),
        downloadButton(outputId = "download_table", label = "Download table")
        ),
        
        #Scatterplots
        column(width = 3,
            #Scatterplot SL
            plotOutput("ScatPlot.SL", hover = hoverOpts(id ="plot_hover_SL"), height = "200px"),
            verbatimTextOutput("hover_info_SL"),
            #Scatterplot 6H
            plotOutput("ScatPlot.6H", hover = hoverOpts(id ="plot_hover_6H"), height = "200px"),
            verbatimTextOutput("hover_info_6H"),
            downloadButton(outputId = "download_plots", label = "Download plots")
            
        ),
        column(width = 3,
            #Scatterplot LB
            plotOutput("ScatPlot.LB", hover = hoverOpts(id ="plot_hover_LB"), height = "200px"),
            verbatimTextOutput("hover_info_LB"),
            #Scatterplot RIL
            plotOutput("ScatPlot.RIL", hover = hoverOpts(id ="plot_hover_RIL"), height = "200px"),
            verbatimTextOutput("hover_info_RIL")
        )
    )
)

# Define server logic
server <- function(input, output) {
        
    output$ScatPlot.SL <- renderPlot({
    Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
    dataset <- input$dataset
    if(dataset==1){
        df <- tmt.SL
        lims <- 1
    }
    else{
        lims <- 2
        df <- cpm.SL
    }
    ScatterplotFun(df=df, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Seedling leaf - Original", lims=lims, Gene=Gene)
    })
    
    output$ScatPlot.RIL <- renderPlot({
        Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
        dataset <- input$dataset
        if(dataset==1){
            df <- tmt.RIL
            lims <- 1
        }
        else{
            lims <- 2
            df <- cpm.RIL
        }
        ScatterplotFun(df=df, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Seedling leaf - RIL hybrids", lims=lims, Gene=Gene)
    })
    
    output$ScatPlot.6H <- renderPlot({
        Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
        dataset <- input$dataset
        if(dataset==1){
            df <- tmt.6H
            lims <- 1
        }
        else{
            lims <- 2
            df <- cpm.6H
        }
        ScatterplotFun(df=df, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Seedling leaf - 6 hybrids", lims=lims, Gene=Gene)
    })
    
    output$ScatPlot.LB <- renderPlot({
        Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
        dataset <- input$dataset
        if(dataset==1){
            df <- tmt.LB
            lims <- 1
        }
        else{
            lims <- 2
            df <- cpm.LB
        }
            ScatterplotFun(df=df, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Leaf blade - Original", lims=lims, Gene=Gene)
    })
    
    output$Cor <- renderTable({
            dataset <- input$dataset
            Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
            if(dataset==1){
            tmt.RIL.6H.Hyb2MP.cor <- read.csv("tmt.RIL.6H.Hyb2MP.cor.csv", stringsAsFactors = F)
            colnames(tmt.RIL.6H.Hyb2MP.cor)[3] <- "Pearson"
            tmt.RIL.6H.Hyb2MP.cor[tmt.RIL.6H.Hyb2MP.cor$Gene %in% Gene,]
            }
            else{
                cpm.RIL.6H.Hyb2MP.cor <- read.csv("cpm.RIL.6H.Hyb2MP.cor.csv", stringsAsFactors = F)
                colnames(cpm.RIL.6H.Hyb2MP.cor)[2] <- "Pearson"
                cpm.RIL.6H.Hyb2MP.cor[cpm.RIL.6H.Hyb2MP.cor$Gene %in% Gene,]
            }
    })
        
    output$hover_info_SL <- renderPrint({
            Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
            dataset <- input$dataset
            tmt.SL.g <- tmt.SL[tmt.SL$Gene %in% Gene,]
            cpm.SL.g <- cpm.SL[cpm.SL$Gene %in% Gene,]
            if(dataset==1){
            if(!is.null(input$plot_hover_SL)){
                hover=input$plot_hover_SL
                dist=sqrt((hover$x-tmt.SL.g$l2.B73xMo17.MP)^2+(hover$y-tmt.SL.g$l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    tmt.SL.g$Accession[which.min(dist)]
                }}}
            else{
                if(!is.null(input$plot_hover_SL)){
                    hover=input$plot_hover_SL
                    dist=sqrt((hover$x-cpm.SL.g$l2.B73xMo17.MP)^2+(hover$y-cpm.SL.g$l2.B73.Mo17)^2)
                    if(min(dist) < .08){
                        cpm.SL.g$Gene[which.min(dist)]
                    }}}
    })
    
    output$hover_info_RIL <- renderPrint({
        Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
        dataset <- input$dataset
        tmt.RIL.g <- tmt.RIL[tmt.RIL$Gene %in% Gene,]
        cpm.RIL.g <- cpm.RIL[cpm.RIL$Gene %in% Gene,]
        if(dataset==1){
            if(!is.null(input$plot_hover_RIL)){
                hover=input$plot_hover_RIL
                dist=sqrt((hover$x-tmt.RIL.g$l2.B73xMo17.MP)^2+(hover$y-tmt.RIL.g$l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    tmt.RIL.g$Accession[which.min(dist)]
                }}}
        else{
            if(!is.null(input$plot_hover_RIL)){
                hover=input$plot_hover_RIL
                dist=sqrt((hover$x-cpm.RIL.g$l2.B73xMo17.MP)^2+(hover$y-cpm.RIL.g$l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    cpm.RIL.g$Gene[which.min(dist)]
                }}}
    })
    
    output$hover_info_6H <- renderPrint({
        Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
        dataset <- input$dataset
        tmt.6H.g <- tmt.6H[tmt.6H$Gene %in% Gene,]
        cpm.6H.g <- cpm.6H[cpm.6H$Gene %in% Gene,]
        if(dataset==1){
            if(!is.null(input$plot_hover_6H)){
                hover=input$plot_hover_6H
                dist=sqrt((hover$x-tmt.6H.g$l2.B73xMo17.MP)^2+(hover$y-tmt.6H.g$l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    tmt.6H.g$Accession[which.min(dist)]
                }}}
        else{
            if(!is.null(input$plot_hover_6H)){
                hover=input$plot_hover_6H
                dist=sqrt((hover$x-cpm.6H.g$l2.B73xMo17.MP)^2+(hover$y-cpm.6H.g$l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    cpm.6H.g$Gene[which.min(dist)]
                }}}
    })
    
    output$hover_info_DN <- renderPrint({
        Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
        tmt.DN.g <- tmt.DN[tmt.DN$Gene %in% Gene,]
        dataset <- input$dataset
        if(dataset==1){
            if(!is.null(input$plot_hover_DN)){
                hover=input$plot_hover_DN
                dist=sqrt((hover$x-tmt.DN.g$day.l2.B73xMo17.MP)^2+(hover$y-tmt.DN.g$day.l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    tmt.DN.g$Accession[which.min(dist)]
                }}}
        else{}
    })
    
    output$hover_info_LB <- renderPrint({
        Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
        dataset <- input$dataset
        tmt.LB.g <- tmt.LB[tmt.LB$Gene %in% Gene,]
        cpm.LB.g <- cpm.LB[cpm.LB$Gene %in% Gene,]
        if(dataset==1){
            if(!is.null(input$plot_hover_LB)){
                hover=input$plot_hover_LB
                dist=sqrt((hover$x-tmt.LB.g$l2.B73xMo17.MP)^2+(hover$y-tmt.LB.g$l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    tmt.LB.g$Accession[which.min(dist)]
                }}}
        else{
            if(!is.null(input$plot_hover_LB)){
                hover=input$plot_hover_LB
                dist=sqrt((hover$x-cpm.LB.g$l2.B73xMo17.MP)^2+(hover$y-cpm.LB.g$l2.B73.Mo17)^2)
                if(min(dist) < .08){
                    cpm.LB.g$Gene[which.min(dist)]
                }}}
    })
    
    output$download_plots <- downloadHandler(
        filename = "ExpressionVis_Plots.pdf",
        content = function(file){
            pdf(file)
            Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
            dataset <- input$dataset
            if(dataset==1){
                df.SL <- tmt.SL
                df.LB <- tmt.LB
                df.6H <- tmt.6H
                df.RIL <- tmt.RIL
                lims <- 1
            }
            else{
                lims <- 2
                df.SL <- cpm.SL
                df.LB <- cpm.LB
                df.6H <- cpm.6H
                df.RIL <- cpm.RIL
            }
            print(ScatterplotFun(df=df.SL, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Seedling leaf - Original", lims=lims, Gene=Gene))
            print(ScatterplotFun(df=df.LB, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Leaf blade - Original", lims=lims, Gene=Gene))
            print(ScatterplotFun(df=df.6H, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Seedling leaf - 6 hybrids", lims=lims, Gene=Gene))
            print(ScatterplotFun(df=df.RIL, xcol="l2.B73xMo17.MP", ycol="l2.B73.Mo17", xlabel="BxM/MP", ylabel="B73/Mo17", title="Seedling leaf - RIL hybrids", lims=lims, Gene=Gene))
            dev.off()
        }
    )
    
    output$download_table <- downloadHandler(
        filename = "ExpresssionVis_PearsonCorTable.csv",
        content = function(file){
            dataset <- input$dataset
            Gene <- as.character(unlist(strsplit(input$Gene,"\\, |\\,| ")))
            if(dataset==1){
                tmt.RIL.6H.Hyb2MP.cor <- read.csv("tmt.RIL.6H.Hyb2MP.cor.csv", stringsAsFactors = F)
                colnames(tmt.RIL.6H.Hyb2MP.cor)[3] <- "Pearson"
                write.csv(tmt.RIL.6H.Hyb2MP.cor[tmt.RIL.6H.Hyb2MP.cor$Gene %in% Gene,], file, row.names = F)
            }
            else{
                cpm.RIL.6H.Hyb2MP.cor <- read.csv("cpm.RIL.6H.Hyb2MP.cor.csv", stringsAsFactors = F)
                colnames(cpm.RIL.6H.Hyb2MP.cor)[2] <- "Pearson"
                write.csv(cpm.RIL.6H.Hyb2MP.cor[cpm.RIL.6H.Hyb2MP.cor$Gene %in% Gene,], file, row.names = F)
            }
        }
    )
}



# Run the application 
shinyApp(ui = ui, server = server)
