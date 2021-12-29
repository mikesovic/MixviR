#' MixviR Shiny Dashboard
#'
#' Open Dashboard To Explore Mutation Data Stored In *samp_mutations* Data Frame Created With `call_mutations()`.
#' @param dates path to optional csv file with cols "SAMP_NAME", "LOCATION", and "DATE". Sample names need to match those in *samp_mutations* data frame created by `call_mutations()`. Dates should be provided in the format *mmddyyyy*.
#' @param lineage.muts path to optional csv file with required cols "Gene", "Mutation", and  "Lineage" containing mutations associated with lineages of interest. See example file at "https://github.com/mikesovic/MixviR/blob/main/mutation_files/outbreak_20211202.csv". Can use this example file by setting lineage.muts to "https://raw.githubusercontent.com/mikesovic/MixviR/main/mutation_files/outbreak_20211202.csv". Two additional columns: "Chr" and "Pos" are optional. These provide the position of the underlying genomic mutation and are used to report the sequencing depths for relevant positions when the mutation of interest is not observed in the sample. Set *get.all.depths* to TRUE. 
#' @param read.muts.from By default, data are read from *samp_mutations* data frame created by `call_mutations()` and written to the global environment. If this data frame has been written to a file (see *write.mut.table* in `call_mutations()`), the dashboard can optionally be opened by providing the path to that file.  
#' @keywords shiny
#' @return Shiny Dashboard to Explore Data
#' @export
#' @examples
#' explore_mutations()

##Which arguments are defined determine what tabs are available in the Shiny dashboard.
##If dates and lineage.def are both defined, get all 4 tabs.
##If lineage.def is defined but dates is not, get plot from dust app and mutation table without dates.
##If dates are defined but lineages aren't, get new mutations table, mutation freq plot, and mutations table (with no lineages).
##If neither are defined, get just mutations table with no lineages.

explore_mutations <- function(dates = NULL, lineage.muts = NULL, read.muts.from = NULL) {
 
  if (is.null(read.muts.from)) {
    samp_data <- samp_mutations %>%
      dplyr::select(SAMP_NAME, CHR, POS, ALT_ID, AF, DP) 
  } else {
    samp_data <- readr::read_tsv(read.muts.from, show_col_types = FALSE) %>%
      dplyr::select(SAMP_NAME, CHR, POS, ALT_ID, AF, DP) 
  }
  
  if (is.null(dates)) {
    
    if (is.null(lineage.muts)) { #get only mutations table with no lineage information column in dashboard
      
      ui <- shiny::fluidPage(shiny::tabsetPanel(shiny::tabPanel("View Mutations",
              shiny::selectInput(inputId = "Sample_mutTable",
                                 label = "Sample",
                                 choices = unique(samp_data$SAMP_NAME)),
                     DT::dataTableOutput("mut_table")
        )
      )
      )
      
      server <- function(input, output, session){
        
        output$mut_table <- DT::renderDataTable({
          
          #create table for selected sample
          samp_data %>% 
            tidyr::separate(col = ALT_ID,
                     into = c("GENE", "MUTATION"),
                     sep = "_") %>%
            dplyr::filter(SAMP_NAME %in% input$Sample_mutTable) %>%
            dplyr::select(SAMP_NAME, CHR, POS, GENE, MUTATION, AF, DP) %>% 
            dplyr::mutate("AF" = round(AF, digits = 3)) %>% 
            dplyr::rename("FREQ" = "AF",
                          "SEQ DEPTH" = "DP") %>%
            as.data.frame()
        })
      }  
    
    ## run Shiny app for no dates/no lineages
    shiny::shinyApp(ui = ui, server = server)
      
      
    } else{ #no dates, but do have lineage defining-mutations - get mutation table and barplot that has all lineages along x-axis.
      
      #read in mutations associated with lineages
      lineage_muts <- readr::read_csv(lineage.muts, col_types = readr::cols_only(Gene = readr::col_character(),
                                                                                 Mutation = readr::col_character(),
                                                                                 Lineage = readr::col_character())) %>%
        tidyr::unite("ALT_ID",
              Gene, Mutation,
              sep = "_") 
      
      #identify mutations that are shared by more than one lineage
      duplicated <- lineage_muts$ALT_ID[duplicated(lineage_muts$ALT_ID)]
      
      #flag mutations that are shared by more than one lineage
      lineage_muts <- lineage_muts %>%
        dplyr::mutate("characteristic" = ifelse(!ALT_ID %in% duplicated, "Y", "N"))
      
      ui <- shiny::fluidPage(shiny::tabsetPanel(
        shiny::tabPanel("Variants Present",
                        shiny::sidebarPanel(width = 2,
                                            shiny::selectInput(inputId = "Sample_VarPres",
                                                              label = "Sample",
                                                              choices = unique(samp_data$SAMP_NAME)),
                                          shiny::sliderInput(inputId = "propThresh",
                                                              label = "Presence Threshold",
                                                              min = 0, max = 1,
                                                              value = 0.5)),
                        shiny::mainPanel(width = 10,
                                         shiny::verticalLayout(
                                           shiny::plotOutput(outputId = "lineages_present"),
                                           shiny::plotOutput(outputId = "lineage_proportions"))
                 )
        ),
        
        shiny::tabPanel("View Mutations",
                 
                 shiny::selectInput(inputId = "Sample_MutTable",
                                   label = "Sample",
                                   choices = unique(samp_data$SAMP_NAME)),
                 DT::dataTableOutput("mut_table")
        )))
      
      server <- function(input, output, session){ 
        
        output$lineages_present <- shiny::renderPlot({
          #get data for the selected Sample
          samp_data <- samp_data %>%
            dplyr::filter(SAMP_NAME %in% input$Sample_VarPres) %>%
            dplyr::filter(AF > 0) #this applies if call_mutations was run with write.all.targets = TRUE - don't want to count things with zero reads.
          
          #get rid of any mutation duplicates, keeping one with highest AF
          #these generally occur when an amino acid change is caused by two or more SNPs - if so, it's repeated for each variant called
          samp_data <- samp_data %>% 
            dplyr::group_by(ALT_ID) %>%
            dplyr::arrange(desc(DP), .by_group = TRUE) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(ALT_ID, .keep_all = TRUE)
          
          #join mutations observed in current sample in to lineage_muts df
          char_muts <- dplyr::left_join(x = lineage_muts, y = samp_data, by = "ALT_ID")
          
          #pull out the mutations that only occur in a single lineage for analysis of lineages present in sample
          char_summary <- char_muts %>%
            dplyr::filter(characteristic == "Y") %>%
            dplyr::group_by(Lineage) %>%
            dplyr::summarize("prop_present" = (sum(!is.na(SAMP_NAME)))/dplyr::n(),
                      "mean_freq" = mean(AF, na.rm = TRUE))
          
          #generate geom_col plot to represent proportion of lineage-characteristic mutations present for each lineage in the sample
          char_summary %>%
            ggplot2::ggplot(ggplot2::aes(x = Lineage, y = prop_present)) +
            ggplot2::geom_col(fill = "red") +
            ggplot2::coord_cartesian(ylim = c(0,1)) +
            ggplot2::theme_classic() +
            ggplot2::ggtitle(paste0("Proportion of Lineage-Characteristic Mutations Present: ", input$Sample_VarPres)) +
            ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.4),
                  axis.text = ggplot2::element_text(size = 13),
                  axis.text.y = ggplot2::element_text(size = 14),
                  plot.title = ggplot2::element_text (color = "black", size= 18, face="bold"),
                  axis.title.y = ggplot2::element_blank()) +
            ggplot2::geom_hline(yintercept = input$propThresh, color = "gray", linetype = 2, alpha = 1)
        })
        
        
        output$lineage_proportions <- shiny::renderPlot({  
          
          #get data for the selected Sample
          samp_data <- samp_data %>%
            dplyr::filter(SAMP_NAME %in% input$Sample_VarPres) %>%
            dplyr::filter(AF > 0) #this applies if call_mutations was run with write.all.targets = TRUE - don't want to count things with zero reads.
          
          #get rid of any mutation duplicates, keeping one with highest AF
          #these generally occur when an amino acid change is caused by two or more SNPs - if so, it's repeated for each variant called
          samp_data <- samp_data %>% 
            dplyr::group_by(ALT_ID) %>%
            dplyr::arrange(desc(DP), .by_group = TRUE) %>%
            dplyr::ungroup() %>%
            dplyr::distinct(ALT_ID, .keep_all = TRUE)
          
          #join mutations observed in current sample in to lineage_muts df
          char_muts <- dplyr::left_join(x = lineage_muts, y = samp_data, by = "ALT_ID")
          
          #pull out the mutations that only occur in a single lineage to estimate relative proportions of each lineage identified as "present" in the sample based on chosen threshold.
          char_summary <- char_muts %>%
            dplyr::filter(characteristic == "Y") %>%
            dplyr::group_by(Lineage) %>%
            dplyr::summarize("prop_present" = (sum(!is.na(SAMP_NAME)))/dplyr::n(),
                      "mean_freq" = mean(AF, na.rm = TRUE)) %>%
            dplyr::mutate_at(dplyr::vars(mean_freq), ~replace(., is.nan(.), 0)) %>%
            dplyr::mutate_at(dplyr::vars(mean_freq), ~replace(., prop_present < input$propThresh, 0))  
          
          #generate geom_col plot to represent proportions of lineages present in the sample.
          #only plotted for lineages identified as "present" based on chosen threshold (otherwise, frequencies were all adjusted to zero above).
          char_summary %>%
            ggplot2::ggplot(ggplot2::aes(x = Lineage, y = mean_freq)) +
            ggplot2::geom_col(fill = "red") +
            ggplot2::coord_cartesian(ylim = c(0,1)) +
            ggplot2::theme_classic() +
            ggplot2::ggtitle(paste0("Estimated Frequency of Variants Present In Sample: ", input$Sample_VarPres)) +
            ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.4),
                  axis.text = ggplot2::element_text(size = 13),
                  plot.title = ggplot2::element_text (color = "black", size= 18, face="bold"),
                  axis.text.y = ggplot2::element_text(size = 14),
                  axis.title.y = ggplot2::element_blank()) +
            ggplot2::geom_hline(yintercept = input$propThresh, color = "gray", linetype = 2, alpha = 1)
          
        })
        
        #generate mutation table that includes column of lineages the mutation is associated with
        output$mut_table <- DT::renderDataTable({
          
          #join lineage muatation info into sample data
          samp_data <- dplyr::left_join(x = samp_data,
                                 y = lineage_muts,
                                 by = "ALT_ID")
          
          #combine all lineages associated with a given mutation - separate with ';'
          #this info will be included as a column in the table
          samp_data <- samp_data %>% 
            dplyr::group_by(SAMP_NAME, CHR, POS, ALT_ID, AF, DP) %>%
            dplyr::summarize("Group" = paste(unique(Lineage),collapse = ";")) %>% 
            dplyr::ungroup()
          
          #create table for selected sample
          samp_data %>% 
            tidyr::separate(col = ALT_ID,
                     into = c("GENE", "MUTATION"),
                     sep = "_") %>%
            dplyr::filter(SAMP_NAME %in% input$Sample_MutTable) %>%
            dplyr::select(SAMP_NAME, CHR, POS, GENE, MUTATION, AF, Group, DP) %>% 
            dplyr::mutate("FREQ" = round(AF, digits = 3)) %>% 
            dplyr::rename("ASSOCIATED LINEAGE(S)" = "Group",
                          "SEQ DEPTH" = "DP") %>%
            as.data.frame()
        })
    }
    
      ## run Shiny app for lineages but no dates
      shiny::shinyApp(ui = ui, server = server)
    }  
    
  } else { #samples have associated dates
    
    #read in file with locations and dates for samples
    dates_df <- readr::read_csv(dates, col_types = "cci")
    
    samp_data <- dplyr::left_join(x = samp_data, y = dates_df, by = "SAMP_NAME") %>%
      dplyr::mutate("date" = lubridate::mdy(DATE)) %>%
      dplyr::mutate("Location" = as.character(LOCATION)) %>%
      dplyr::select(-DATE, -LOCATION)
    
    if (is.null(lineage.muts)) { #have dates but no lineage definitions - output new mutations table, mutation freq plot, and mutations table (without associated lineages column).
      
      ui <- shiny::fluidPage(shiny::tabsetPanel(
        shiny::tabPanel("New Mutations",
                        shiny::selectInput(inputId = "MaxDate",
                             label = "Mutations First Detected On Or After...",
                             choices = unique(samp_data$date)),
                        shiny::numericInput(inputId = "DPThresh",
                              label = "Minimum Depth",
                              value = 5),
                 DT::dataTableOutput("new_muts")),
        shiny::tabPanel("Mutation Frequencies",
                        shiny::textInput(inputId = "Muts",
                           label = "Mutation(s) (i.e S_D614G; separate with commas)"),
                        shiny::selectInput(inputId = "Location_Muts",
                             label = "Location",
                             choices = unique(samp_data$Location),
                             multiple = TRUE),
                        shiny::plotOutput(outputId = "mut_freqs")),
        shiny::tabPanel("View Mutations",
                       shiny::selectInput(inputId = "Sample_MutTable",
                             label = "Location",
                             choices = unique(samp_data$Location)),
                        shiny::selectInput(inputId = "Date",
                             label = "Date",
                             choices = unique(samp_data$date)),
                 DT::dataTableOutput("mut_table")
        )
      ))
      
      server <- function(input, output, session){
        
        #create table with mutations first observed on or after a selected date
        output$new_muts <- DT::renderDataTable({
          
          #filter for mutations observed with a specific sequencing depth
          samp_data <- samp_data %>%
            dplyr::filter(DP >= input$DPThresh) %>%
            dplyr::filter(AF > 0) #this applies if call_mutations was run with write.all.targets = TRUE - don't want to count things with zero reads.
          
          #ID the mutations observed prior to the chosen date - these will be filtered out
          already_observed <- samp_data %>%
            dplyr::filter(date < input$MaxDate) %>%
            dplyr::select(ALT_ID, AF, SAMP_NAME) %>%
            dplyr::distinct(ALT_ID, .keep_all = TRUE) %>%
            unlist()
          
          #remove duplicate mutations (occur if multiple underlying mutations give rise to specific amino acid change)
          all_uniqs <- samp_data %>%
            dplyr::select(ALT_ID, AF, SAMP_NAME) %>%
            dplyr::distinct(ALT_ID, .keep_all = TRUE) %>%
            unlist()
          
          #ID mutations first observed on or after selected date
          new_mutations <- all_uniqs[!all_uniqs %in% already_observed]
          
          #filter for new mutations
          samp_data %>% dplyr::filter(ALT_ID %in% new_mutations)
        })
        
        output$mut_freqs <- shiny::renderPlot({
          #get set of mutations of interest
          #these are entered as comma separated values in text box. Allows for spaces between or not.
          target_muts <- stringr::str_split(input$Muts, pattern = ",\\s?") %>% unlist()
          
          #plot frequencies of selected mutations for the selected Sample(s)
          #multiple mutations are differentiated by color
          #multiple sites are plotted as separate facets
          samp_data %>%
            dplyr::filter(Location %in% input$Location_Muts) %>%
            dplyr::filter(ALT_ID %in% target_muts) %>%
            ggplot2::ggplot(ggplot2::aes(x = date, y = AF, color = Location)) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::facet_wrap(dplyr::vars(c(ALT_ID))) +
            ggplot2::coord_cartesian(ylim = c(0,1)) +
            ggplot2::theme(legend.text = ggplot2::element_text(size = 13),
                           legend.title = ggplot2::element_text(size = 14),
                           strip.text = ggplot2::element_text(size = 14),
                           axis.text.x = ggplot2::element_text(size = 12),
                           axis.text.y = ggplot2::element_text(size = 13),
                           axis.title.y = ggplot2::element_text(size = 13)) +
            ggplot2::labs(y = "Mutation Frequency", x = NULL)
        })
        
        #update the "Date" drop-down box to only include dates for which samples exist for the selected Sample
        shiny::observeEvent(input$Sample_MutTable, {
          date_options <- samp_data %>%
            dplyr::filter(Location == input$Sample_MutTable) %>%
            dplyr::distinct(date) 
          date_options <- date_options$date
          
          shiny::updateSelectInput(session, "Date",
                            label = "Date",
                            choices = date_options)
        })
        
        #create table with all mutations observed in selected sample. Will not include column showing any lineages they have been associated with.
        output$mut_table <- DT::renderDataTable({
          
          #create table for selected sample/date
          samp_data %>% 
            tidyr::separate(col = ALT_ID,
                     into = c("GENE", "MUTATION"),
                     sep = "_") %>%
            dplyr::filter(Location %in% input$Sample_MutTable) %>%
            dplyr::filter(date == input$Date) %>%
            dplyr::select(SAMP_NAME, Location, date, CHR, POS, GENE, MUTATION, AF) %>% 
            dplyr::mutate("FREQ" = round(AF, digits = 3)) %>% 
            dplyr::rename("DATE" = "date",
                          "LOCATION" = "Location") %>%
            dplyr::select(-AF) %>%
            dplyr::select(SAMP_NAME, LOCATION, DATE, CHR, POS, GENE, MUTATION, FREQ) %>%
            as.data.frame()
        })
        
      }  
      
      ## run Shiny app with dates but no lineages
      shiny::shinyApp(ui = ui, server = server)
      
    } else { #have both dates and lineage definitions - output all tabs
      
      #get mutations associated with lineages
      lineage_muts <- readr::read_csv(lineage.muts, col_types = readr::cols_only(Gene = readr::col_character(),
                                                                                 Mutation = readr::col_character(),
                                                                                 Lineage = readr::col_character())) %>%
        tidyr::unite("ALT_ID",
                     Gene, Mutation,
                     sep = "_") 
      
      #identify mutations that are shared by more than one lineage
      duplicated <- lineage_muts$ALT_ID[duplicated(lineage_muts$ALT_ID)]
      
      #flag mutations that are shared by more than one lineage
      lineage_muts <- lineage_muts %>%
        dplyr::mutate("characteristic" = ifelse(!ALT_ID %in% duplicated, "Y", "N"))
      
      ui <- shiny::fluidPage(shiny::tabsetPanel(
        shiny::tabPanel("Variants Present",
                        shiny::sidebarPanel(width = 2,
                                            shiny::selectInput(inputId = "Location",
                                              label = "Location",
                                              choices = unique(samp_data$Location)),
                                            shiny::checkboxGroupInput(inputId = "Lineages",
                                                 label = "Variant/Lineage",
                                                 choices = unique(lineage_muts$Lineage)),
                                            shiny::selectInput(inputId = "DepthTag",
                                                  label = "Calculate Mean Depths On...",
                                                  choices = c("All Mutations", "Lineage-Characteristic Mutations Only")),
                                            shiny::sliderInput(inputId = "propThresh",
                                                label = "Presence Threshold",
                                                min = 0, max = 1,
                                                value = 0.5),
                                            shiny::selectInput(inputId = "PlotType",
                                                               label = "View Plot As...",
                                                               choices = c("Line Plot", "Area Plot"))),
                        shiny::mainPanel(width = 10,
                                         shiny::verticalLayout(
                                           shiny::plotOutput(outputId = "lineages_present"),
                                           shiny::plotOutput(outputId = "lineage_proportions"))
                        )
        ),
        
        shiny::tabPanel("New Mutations",
                        shiny::selectInput(inputId = "MaxDate",
                             label = "Mutations First Detected On Or After...",
                             choices = unique(samp_data$date)),
                        shiny::numericInput(inputId = "DPThresh",
                              label = "Minimum Depth",
                              value = 5),
                 DT::dataTableOutput("new_muts")
                 
        ),
        shiny::tabPanel("Mutation Frequencies",
                        shiny::textInput(inputId = "Muts",
                           label = "Mutation(s) (i.e S_D614G; separate with commas)"),
                        shiny::selectInput(inputId = "Location_Muts",
                             label = "Location",
                             choices = unique(samp_data$Location),
                             multiple = TRUE),
                        shiny::plotOutput(outputId = "mut_freqs")
                 
        ),
        
        shiny::tabPanel("View Mutations",
                 
                        shiny::selectInput(inputId = "Location_MutTable",
                             label = "Location",
                             choices = unique(samp_data$Location)),
                        shiny::selectInput(inputId = "Date",
                             label = "Date",
                             choices = unique(samp_data$date)),
                 DT::dataTableOutput("mut_table")
        )
      ))
      
      server <- function(input, output, session){
        
        output$lineages_present <- shiny::renderPlot({
          #filter for data for the selected Location
          samp_data <- samp_data %>%
            dplyr::filter(Location %in% input$Location) %>%
            dplyr::filter(AF > 0) #this applies if call_mutations was run with write.all.targets = TRUE - don't want to count things with zero reads.
          
          #create master df that will store data for all selected lineages/variants
          all_summary <- data.frame()
          
          #loop over each selected lineage
          for (i in 1:length(input$Lineages)){
            #get the mutations characteristic of the current lineage
            muts <- lineage_muts %>%
              dplyr::filter(characteristic == "Y") %>%
              dplyr::filter(Lineage %in% input$Lineages[[i]]) %>%
              dplyr::select(ALT_ID) %>%
              unlist()
            
            #get the number of mutations characteristic of the current lineage
            lineage_n <- length(muts)
            
            #get proportion of characteristic mutations observed on each date for current Location
            summary <- samp_data %>%
              dplyr::filter(ALT_ID %in% muts) %>%
              #remove any duplicated mutations, keeping one with highest read depth
              dplyr::group_by(SAMP_NAME, date, ALT_ID) %>%
              dplyr::arrange(desc(DP), .by_group = TRUE) %>%
              dplyr::ungroup() %>%
              dplyr::distinct(SAMP_NAME, date, ALT_ID, .keep_all = TRUE) %>%
              dplyr::group_by(SAMP_NAME, date) %>%
              dplyr::summarise("prop_observed" = dplyr::n()/lineage_n,
                        "lineage" = input$Lineages[[i]],
                        "lineage_n" = paste0(input$Lineages[[i]], " (", lineage_n, ")"),
                        "Mean Coverage" = mean(DP)) %>%
              dplyr::ungroup() %>% 
              dplyr::mutate("ID" = paste0(SAMP_NAME, "_", date))
            
            #recalculate Mean Coverage based on all mutations and add in to 'summary' from above
            summary2 <- samp_data %>%
              dplyr::group_by(SAMP_NAME, date, ALT_ID) %>%
              dplyr::arrange(desc(DP), .by_group = TRUE) %>%
              dplyr::ungroup() %>%
              dplyr::distinct(SAMP_NAME, date, ALT_ID, .keep_all = TRUE) %>%
              dplyr::group_by(SAMP_NAME, date) %>%
              dplyr::summarise("Mean_Coverage_All" = mean(DP)) %>%
              dplyr::ungroup() %>%
              dplyr::mutate("ID" = paste0(SAMP_NAME, "_", date)) %>%
              dplyr::select(ID, Mean_Coverage_All)
            
            summary <- dplyr::left_join(x = summary, y = summary2, by = "ID")
            
            #add data for current lineage to master df
            all_summary <- dplyr::bind_rows(all_summary, summary)
            
          }  
          
          #use days since earliest sample date as x axis (lubridate interval)
          ints <- lubridate::interval(start = min(all_summary$date), end = all_summary$date)
          
          #convert seconds to days
          all_summary <- all_summary %>%
            dplyr::mutate("days" = lubridate::int_length(ints)/86400)
          
          #get vectors for labeling x-axis
          break_dates <- unique(all_summary$date)
          break_days <- unique(all_summary$days)
          
          #check which type of mean depth to plot
          if (input$DepthTag == "All Mutations") {
            
            #replace Mean Coverage col with data from Mean_Coverage_All
            all_summary <- all_summary %>% 
              dplyr::select(-`Mean Coverage`) %>%
              dplyr::rename("Mean Coverage" = "Mean_Coverage_All")
          }
          
          #generate plot
          if (input$PlotType == "Area Plot") {
            all_summary %>% 
              dplyr::rename("Lineage (n)" = "lineage_n") %>%
              ggplot2::ggplot(ggplot2::aes(x = days, y = prop_observed)) +
              #ggplot2::geom_line(ggplot2::aes(colour = `Lineage (n)`, group = `Lineage (n)`)) +
              ggplot2::geom_area(ggplot2::aes(fill = `Lineage (n)`), stat = "identity", position = "identity", alpha = 0.5) +
              ggplot2::coord_cartesian(ylim = c(0,1)) +
              ggplot2::theme_classic() +
              ggplot2::ggtitle(paste0("Proportion of Lineage-Characteristic Mutations Present: ", input$Location)) +
              ggplot2::scale_x_continuous(breaks = break_days, labels = break_dates) +
              ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                    axis.text.x = ggplot2::element_text(angle = 45, size = 12, vjust = 0.4),
                    axis.text.y = ggplot2::element_text(size = 13),
                    plot.title = ggplot2::element_text (color = "black", size= 15, face="bold"),
                    axis.title.y = ggplot2::element_blank(),
                    legend.text = ggplot2::element_text(size = 13),
                    legend.title = ggplot2::element_text(size = 14)) +
              ggplot2::geom_hline(yintercept = input$propThresh, color = "red", linetype = 2, alpha = 0.6) +
              ggplot2::geom_point(data = all_summary, mapping = ggplot2::aes(x = days, 
                                                           y = prop_observed,
                                                           color = `Mean Coverage`), size = 2) +
              ggplot2::scale_color_gradient(low = "#56B1F7", 
                                    high = "#132B43",
                                    limits = c(0,4200))
          } else {
            all_summary %>% 
              dplyr::rename("Lineage (n)" = "lineage_n") %>%
              ggplot2::ggplot(ggplot2::aes(x = days, y = prop_observed)) +
              ggplot2::geom_line(ggplot2::aes(colour = `Lineage (n)`, group = `Lineage (n)`)) +
              ggplot2::coord_cartesian(ylim = c(0,1)) +
              ggplot2::theme_classic() +
              ggplot2::ggtitle(paste0("Proportion of Lineage-Characteristic Mutations Present: ", input$Location)) +
              ggplot2::scale_x_continuous(breaks = break_days, labels = break_dates) +
              ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                             axis.text.x = ggplot2::element_text(angle = 45, size = 12, vjust = 0.4),
                             axis.text.y = ggplot2::element_text(size = 13),
                             plot.title = ggplot2::element_text (color = "black", size= 15, face="bold"),
                             axis.title.y = ggplot2::element_blank(),
                             legend.text = ggplot2::element_text(size = 13),
                             legend.title = ggplot2::element_text(size = 14)) +
              ggplot2::geom_hline(yintercept = input$propThresh, color = "red", linetype = 2, alpha = 0.6) +
              ggplot2::geom_point(data = all_summary, mapping = ggplot2::aes(x = days, 
                                                                             y = prop_observed,
                                                                             fill = `Mean Coverage`), size = 2) +
              ggplot2::scale_fill_gradient(low = "#56B1F7", 
                                            high = "#132B43",
                                            limits = c(0,4200))
            
          }
          
          
        })
        
        
        output$lineage_proportions <- shiny::renderPlot({
          #filter for data for the selected Location
          samp_data <- samp_data %>%
            dplyr::filter(Location %in% input$Location) %>%
            dplyr::filter(AF > 0) #this applies if call_mutations was run with write.all.targets = TRUE - don't want to count things with zero reads.
          
          #create master df that will store data for all selected lineages
          all_summary <- data.frame()
          
          #loop over each selected lineage
          for (i in 1:length(input$Lineages)){
            
            #get the mutations characteristic of the current lineage
            muts <- lineage_muts %>%
              dplyr::filter(characteristic == "Y") %>%
              dplyr::filter(Lineage %in% input$Lineages[[i]]) %>%
              dplyr::select(ALT_ID) %>%
              unlist()
            
            #get the number of mutations characteristic of the current lineage
            lineage_n <- length(muts)
            
            #get proportion of characteristic mutations observed in current sample and avg frequency of observed mutations
            summary <- samp_data %>%
              dplyr::filter(ALT_ID %in% muts) %>%
              dplyr::group_by(SAMP_NAME, date, ALT_ID) %>%
              dplyr::arrange(desc(DP), .by_group = TRUE) %>%
              dplyr::ungroup() %>%
              dplyr::distinct(SAMP_NAME, date, ALT_ID, .keep_all = TRUE) %>%
              dplyr::group_by(SAMP_NAME, date) %>%
              dplyr::summarise("prop_observed" = dplyr::n()/lineage_n,
                        "lineage" = input$Lineages[[i]],
                        "lineage_n" = paste0(input$Lineages[[i]], " (", lineage_n, ")"),
                        "Mean Coverage" = mean(DP),
                        "avg_freq" = mean(AF))
            
            #add data for current lineage to master df
            all_summary <- dplyr::bind_rows(all_summary, summary)
            
          }  
          
          #use days since earliest sample date as x axis (lubridate interval)
          ints <- lubridate::interval(start = min(all_summary$date), end = all_summary$date)
          
          #convert seconds to days
          all_summary <- all_summary %>%
            dplyr::mutate("days" = lubridate::int_length(ints)/86400)
          
          #get vectors for labeling x-axis
          break_dates <- unique(all_summary$date)
          break_days <- unique(all_summary$days)
          
          #generate plot
          
          #only plot bars for lineages defined as "present" in the sample based on chosen proportion threshold 
          all_summary$avg_freq[all_summary$prop_observed < input$propThresh] <- 0
          
          all_summary %>%
            dplyr::rename("Lineage (n)" = "lineage_n") %>%
            ggplot2::ggplot(ggplot2::aes(x = days, y = avg_freq, fill = `Lineage (n)`)) +
            ggplot2::geom_bar(stat = "identity",
                     position = "stack") +
            ggplot2::theme_classic() +
            ggplot2::ggtitle(paste0("Estimated Frequency of Each Lineage Present: ", input$Location)) +
            ggplot2::scale_x_continuous(breaks = break_days, labels = break_dates) +
            ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                  axis.text.x = ggplot2::element_text(angle = 45, vjust = 0.4),
                  axis.text = ggplot2::element_text(size = 13),
                  plot.title = ggplot2::element_text (color = "black", size= 15, face="bold"),
                  axis.title.y = ggplot2::element_blank(),
                  legend.text = ggplot2::element_text(size = 13),
                  legend.title = ggplot2::element_text(size = 14))
          
        })
        
        #create table with mutations first observed on or after a selected date
        output$new_muts <- DT::renderDataTable({
          
          #merge in lineage information
          samp_data <- dplyr::left_join(x = samp_data,
                                        y = lineage_muts,
                                        by = "ALT_ID") %>%
            dplyr::filter(AF > 0) #this applies if call_mutations was run with write.all.targets = TRUE - don't want to count things with zero reads.
          
          #filter for mutations observed with a specific sequencing depth
          #get all lineages a given mutation is associated with as a ';'-delimited string - will be printed to column
          samp_data <- samp_data %>%
            dplyr::filter(DP >= input$DPThresh) %>%
            dplyr::group_by(SAMP_NAME, date, ALT_ID, AF) %>%
            dplyr::summarize("Group" = paste(unique(Lineage),collapse = ";")) %>% 
            dplyr::ungroup() 
    
          #ID the mutations observed prior to the chosen date - these will be filtered out
          already_observed <- samp_data %>%
            dplyr::filter(date < input$MaxDate) %>%
            dplyr::select(ALT_ID, AF, SAMP_NAME, Group) %>%
            dplyr::distinct(ALT_ID, .keep_all = TRUE) %>%
            unlist()
          
          #remove duplicate mutations (occur if multiple underlying mutations give rise to specific amino acid change)
          all_uniqs <- samp_data %>%
            dplyr::select(ALT_ID, AF, SAMP_NAME, Group) %>%
            dplyr::distinct(ALT_ID, .keep_all = TRUE) %>%
            unlist()
          
          #ID mutations first observed on or after selected date
          new_mutations <- all_uniqs[!all_uniqs %in% already_observed]
          
          samp_data %>% dplyr::filter(ALT_ID %in% new_mutations)
        })
        
        
        output$mut_freqs <- shiny::renderPlot({
          #get set of mutations of interest
          #these are entered as comma separated values in text box. Allows for spaces between or not.
          target_muts <- stringr::str_split(input$Muts, pattern = ",\\s?") %>% unlist()
          
          #plot frequencies of selected mutations for the selected Sample(s)
          #multiple mutations are differentiated by color
          #multiple sites are plotted as separate facets
          samp_data %>%
            dplyr::filter(Location %in% input$Location_Muts) %>%
            dplyr::filter(ALT_ID %in% target_muts) %>%
            ggplot2::ggplot(ggplot2::aes(x = date, y = AF, color = Location)) +
            ggplot2::geom_line() +
            ggplot2::geom_point(size = 3) +
            ggplot2::facet_wrap(dplyr::vars(c(ALT_ID))) +
            ggplot2::coord_cartesian(ylim = c(0,1)) +
            ggplot2::theme(legend.text = ggplot2::element_text(size = 13),
                           legend.title = ggplot2::element_text(size = 14),
                           strip.text = ggplot2::element_text(size = 14),
                           axis.text.x = ggplot2::element_text(size = 12),
                           axis.text.y = ggplot2::element_text(size = 13),
                           axis.title.y = ggplot2::element_text(size = 13)) +
            ggplot2::labs(y = "Mutation Frequency", x = NULL)
        })
        
        #update the "Date" drop-down box to only incude dates for which samples exist for the selected Sample
        shiny::observeEvent(input$Location_MutTable, {
          date_options <- samp_data %>%
            dplyr::filter(Location == input$Location_MutTable) %>%
            dplyr::distinct(date) 
          date_options <- date_options$date
          
          shiny::updateSelectInput(session, "Date",
                            label = "Date",
                            choices = date_options)
        })
        
        #create table with all mutations observed in selected sample and any lineages they have been associated with.
        output$mut_table <- DT::renderDataTable({
          
          #merge in info on lineages associated with specific mutations
          samp_data <- dplyr::left_join(x = samp_data,
                                 y = lineage_muts,
                                 by = "ALT_ID")
          
          #combine all lineages associated with a given mutation - separate with ';'
          samp_data <- samp_data %>% 
            dplyr::group_by(SAMP_NAME, Location, date, CHR, POS, ALT_ID, AF, DP) %>%
            dplyr::summarize("Group" = paste(unique(Lineage),collapse = ";")) %>%
            dplyr::ungroup()
          
          #create table for selected sample/date
          samp_data %>% 
            tidyr::separate(col = ALT_ID,
                     into = c("GENE", "MUTATION"),
                     sep = "_") %>%
            dplyr::filter(Location %in% input$Location_MutTable) %>%
            dplyr::filter(date == input$Date) %>%
            dplyr::select(SAMP_NAME, Location, date, CHR, POS, GENE, MUTATION, AF, Group, DP) %>% 
            dplyr::mutate("AF" = round(AF, digits = 3)) %>% 
            dplyr::rename("DATE" = "date", 
                          "ASSOCIATED LINEAGE(S)" = "Group",
                          "FREQ" = "AF",
                          "SEQ DEPTH" = "DP",
                          "LOCATION" = "Location") %>%
            as.data.frame()
        })
        
      } 
      ## run Shiny app for lineages and dates
      shiny::shinyApp(ui = ui, server = server)
      
    }
  }
}



