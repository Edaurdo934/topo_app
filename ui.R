library(shiny)
library(shinyauthr)
library(shinythemes)
library(dplyr)
library(lubridate)
library(RPostgreSQL)
library(DBI)
library(pool)
library(DT)
library(shinyjs)

# Días para que la sesión expire
cookie_caducidad<- 5

# UI
shinyUI(
    
    bootstrapPage(
        # UI panel
        shinyjs::useShinyjs(),
        tags$head(
            ## Favicon
            tags$head(
                tags$link(rel="shortcut icon", href="favicon-16x16.png")
            ),
            # Note the wrapping of the string in HTML()
            tags$style(HTML(
                "
                #modelo, #texto_corregido {font-size:12px; font-style:italic; 
                           overflow-y:scroll; max-height: 500px; background: ghostwhite;
                
                "
            ))
        ),
        shinyauthr::loginUI(id = "login", 
                            cookie_expiry = cookie_caducidad,
                            title= div(img(src="SG Logo.png", align="center", style="width:60%; height:60%;"),
                                       h2("Welcome")),
                            additional_ui = fluidRow(
                                column(10,
                                       actionLink("registro",label = "Register", class="pull-right"),
                                       a("Help", href = " mailto:soporte@solucionesgeograficas.pe", target = "_blank")
                                       )
                            )
        ),
        theme = shinytheme("cerulean"),
        uiOutput("navbar"),
        # Logout button, cambiar esto y ponerlo en un luagar mejor
        div(shinyauthr::logoutUI(id = "logout", label="Log Out"), style="position: absolute; top: 8px; right: 80px; z-index: 1000;"),
        
        HTML('<footer style="background-color:#e3e3e3; height:50px; width:100%;margin-top:100px;">
                      <p>&copy; 2022 Soluciones Geográficas - all rights reserved</p>
                    </footer>')
        
    )
)
