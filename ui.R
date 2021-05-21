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

# UI
shinyUI(
    
    bootstrapPage(
        # UI panel
        shinyjs::useShinyjs(),
        tags$head(
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
                            title= "Bienvenido",
                            error_message = "Usuario o contraseña incorrectas",
                            additional_ui = fluidRow(
                                column(10,
                                       actionLink("registro",label = "Registrate", class="pull-right"),
                                       a("Ayuda", href = " mailto:jose.oro.joj@gmail.com", target = "_blank")
                                       )
                            )
        ),
        theme = shinytheme("cerulean"),
        uiOutput("navbar"),
        # Logout button, cambiar esto y ponerlo en un luagar mejor
        div(shinyauthr::logoutUI(id = "logout", label="Salir"), style="position: absolute; top: 8px; right: 80px; z-index: 1000;")
        
    )
)
