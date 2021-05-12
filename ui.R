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
        shinyauthr::loginUI(id = "login", 
                            cookie_expiry = cookie_caducidad,
                            title= "Bienvenido",
                            error_message = "Usuario o contraseña incorrectas",
                            additional_ui = fluidRow(
                                column(10,
                                       actionLink("registro",label = "Registrate", class="pull-right"),
                                       a("Ayuda", href = " mailto:hector.br934@gmail.com", target = "_blank")
                                       )
                            )
        ),
        theme = shinytheme("cerulean"),
        uiOutput("navbar"),
        # Logout button, cambiar esto y ponerlo en un luagar mejor
        div(shinyauthr::logoutUI(id = "logout"), style="position: absolute; top: 8px; right: 80px; z-index: 1000;")
        
    )
)
