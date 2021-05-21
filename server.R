library(shiny)
library(shinyjs)
library(shinyauthr)
library(dplyr)
library(plyr)
library(sodium)
library(lubridate)
library(RPostgreSQL)
library(DBI)
library(pool)
library(sf)
library(leaflet)
### Comentario
#La tabla user_table muestra la información de cada usuario

## Función que calcula azimuth entre puntos
st_azimuth = function(x, y) {
    
    # Checks
    stopifnot(all(st_is(x, "POINT")))
    stopifnot(all(st_is(y, "POINT")))
    
    # Extract geometry
    x = st_geometry(x)
    y = st_geometry(y)
    
    # Recycle 'x' or 'y' if necessary
    if(length(x) < length(y)) {
        ids = rep(1:length(x), length.out = length(y))
        x = x[ids]
    }
    if(length(y) < length(x)) {
        ids = rep(1:length(y), length.out = length(x))
        y = y[ids]
    }
    
    # Get coordinate matrices
    x_coords = st_coordinates(x)
    y_coords = st_coordinates(y)
    
    # Calculate azimuths
    x1 = x_coords[, 1]
    y1 = x_coords[, 2]
    x2 = y_coords[, 1]
    y2 = y_coords[, 2]
    az = (180 / pi) * atan2(x2 - x1, y2 - y1)
    names(az) = NULL
    az[az < 0] = az[az < 0] + 360
    
    # Replace with 'NA' for identical points
    az[x1 == x2 & y1 == y2] = NA
    
    # Return
    return(az)
    
}

options(digits=12, scipen=999)
# Funciones matriciales
## Cálculo matricial Ajuste local horizontal
###Función que recibe una matriz con columnas X,Y,E,N y devuelve las matrices A, X, L y V
crear_X <- function(x){
    matriz_a <- matrix(nrow = 2*nrow(x), ncol = 4)
    matriz_L <- matrix(nrow = 2*nrow(x), ncol = 1)
    
    for(i in 1: nrow(x)){
        matriz_a[(2*i-1), 1] <- x[i, 1]
        matriz_a[(2*i-1),2] <- -x[i,2]
        matriz_a[(2*i-1), 3] <- 1
        matriz_a[(2*i-1), 4] <- 0
        matriz_a[(2*i), 1] <- x[i, 2]
        matriz_a[(2*i), 2] <- x[i, 1]
        matriz_a[(2*i), 3] <- 0
        matriz_a[(2*i), 4]<- 1
    }
    for (i in 1:nrow(x)) {
        matriz_L[(2*i-1),1] <- x[i,3]
        matriz_L[(2*i),1] <- x[i,4]
        
    }
    matriz_v <-matriz_a%*%(solve(t(matriz_a) %*% diag(nrow(matriz_a)) %*% matriz_a, tol = 1e-25) %*% t(matriz_a) %*% diag(nrow(matriz_a)) %*% matriz_L)-matriz_L
    matriz_x <- solve(t(matriz_a) %*% diag(nrow(matriz_a)) %*% matriz_a, tol = 1e-25) %*% t(matriz_a) %*% diag(nrow(matriz_a)) %*% matriz_L
    
    matrices<-list(matriz_a,matriz_x ,matriz_L, matriz_v)
    return(matrices)
    
    
}
###Función que hace los ajustes y arroja los nuevos valores para el cálculo local horizontal. Recibe una matriz X Y  y regresa una matriz E N
resultados <- function(matriz_elementos, x){
    a <- x[1, 1]
    b <- x[2, 1]
    T_x <- x[3, 1]
    T_y <- x[4, 1]
    
    p <- matrix(nrow = nrow(matriz_elementos), ncol = 2)
    
    for (i in 1:nrow(matriz_elementos)) {
        p[i,1] <- a*matriz_elementos[i,1]-b*matriz_elementos[i,2]+T_x
        p[i,2] <- b*matriz_elementos[i,1]+a*matriz_elementos[i,2]+T_y
        
    }
    return(p)
    
}

## Cálculo matricial ajuste local vertical
### Recibe una matriz E,N, H, h y devuelve las matrices X el vector Y y el vetor B
crear_B <- function(x,h){
    matriz_x <- matrix(nrow = nrow(x), ncol = 3)
    for (i in 1:nrow(x)) {
        matriz_x[i,1] <- 1
        matriz_x[i,2] <- x[i,1]
        matriz_x[i,3] <- x[i,2]
    }
    
    y <- matrix(nrow = nrow(x), ncol = 1)
    
    for (j in 1:nrow(x)) {
        y[j,1] <- x[j,4]-x[j,3] 
    }
    
    matriz_trans_x <- t(matriz_x)
    mult_tx_x <- matriz_trans_x %*% matriz_x
    mat_inv <- solve(mult_tx_x, tol = 1e-25)
    tx_x_y <- matriz_trans_x %*% y
    
    mat_b <- mat_inv %*% tx_x_y
    
    matrices<-list(matriz_x,y,mat_b)
    
    return(matrices)
    
}
### Recibe la matriz E,N,h y devuelve el vector ondulación geoidal y H 
resultados_H <- function(x, mat_b){
    ondulacion <- matrix(nrow = nrow(x), ncol = 1)
    for (i in 1:nrow(x)) {
        ondulacion[i, 1] <- mat_b[1,1]+mat_b[2,1]*x[i,1]+mat_b[3,1]*x[i,2]
        
    }
    
    h <- matrix(nrow = nrow(x), ncol = 1)
    
    for (i in 1:nrow(x)) {
        h[i,1] <- x[i,3]
        
    }
    
    matriz_h <- matrix(nrow = nrow(x), ncol = 1)
    
    for (i in 1:nrow(x)) {
        matriz_h[i,1] <- h[i,1]-ondulacion[i,1]
    }
    
    resultados_f<-list(ondulacion, matriz_h)
    
    return(resultados_f)
}


##### función para el ajuste UTM_plano
## Recibe un data.frame con columnas E,N y h  además del CRS devuelve un data.frame con "Lat" "Long"  "k_esc" "K_ele" 
#"K_com"       "dist_UTM"    "Acimut_rad"  "acimut_grad" "dist_top"    "E_top"       "N_top" 
funcion_UTM_planas <- function(p, crs_input){
    a <- 6378137.00
    b <- 6356752.314
    e2<-0.00669438
    et2 <- 0.006739497
    c <- 6399593.626
    k_0 <- 0.9996
    
    h <- st_as_sf(p, coords = c("E","N"), crs=crs_input)
    q <- st_transform(h, crs = 4326)
    m <- st_coordinates(q)
    
    
    k_esc <- matrix(nrow = nrow(m), ncol = 1)
    k_ele <- matrix(nrow = nrow(m), ncol = 1)
    k_com <- matrix(nrow = nrow(m), ncol = 1)
    azimut_grad <- matrix(nrow = nrow(m), ncol = 1)
    azimut_rad <- matrix(nrow = nrow(m), ncol = 1)
    distance_UTM <- matrix(nrow = nrow(m), ncol = 1)
    distance_top <- matrix(nrow = nrow(m), ncol = 1)
    e_top <- matrix(nrow = nrow(m), ncol = 1)
    n_top <- matrix(nrow = nrow(m), ncol = 1)
    
    for (i in 1:nrow(m)) {
        
        x <- 500000-p[i,1]
        q <- 0.000001*x
        
        n_k <- a/((1-e2*(sin(m[i,2]*pi/180))^2)^(1/2))
        p_k <- (10^(12))*((1+et2*(cos(m[i,2]*pi/180))^2)/(2*n_k^2*k_0^2))
        k_esc[i,1] <- k_0*(1+p_k*q^2+0.00003*q^4)
        k_ele[i,1] <- (a*(1-e2))/(a*(1-e2)+p[i,3]*(1-e2*(sin(m[i,2]*pi/180))^2)^(3/2))
        k_com[i,1] <- k_esc[i,1]*k_ele[i,1]
        if(i==1){
            azimut_grad[i,1] <- 90
        }
        else{
            azimut_grad[i,1] <- st_azimuth(h[1,1], h[i,1])
        }
        azimut_rad[i,1] <- azimut_grad[i,1]*pi/180
        distance_UTM[i,1] <- st_distance(h)[i,1]
        distance_top[i,1] <- distance_UTM[i,1]/mean(c(k_com[i,1], k_com[1,1]))
        e_top[i,1] <- p[1,1]+distance_top[i,1]*sin(azimut_rad[i,1])
        n_top[i,1] <- p[1,2]+distance_top[i,1]*cos(azimut_rad[i,1])
    }
    
    datos <- data.frame("Lat"=m[,2], "Long"=m[,1], "k_esc"=k_esc, "K_ele"=k_ele, "K_com"=k_com, "dist_UTM"=distance_UTM, "Acimut_rad"=azimut_rad, "acimut_grad"=azimut_rad, "dist_top"=distance_top,"E_top"=e_top, "N_top"=n_top)
    return(datos)
    
    
    
}

############### Conexión a la base de datos de aws y ajuste de parámetros del login########
conexion_base <- dbPool(
    drv = dbDriver("PostgreSQL", max.con = 100),
    dbname = "topo_app",
    host = "database-topoapp.cic5fcsqtcdl.us-east-2.rds.amazonaws.com",
    user = "postgres",
    password = "postgres_jose.oro",
    idleTimeout = 3600000
)



######################## Servidor################################
shinyServer(function(input, output, session) {
    ## Arregla el bug del cambio de vista al momento de cargar datos para corrección proceso local
    onclick("archivo_corregir",
            runjs("
                  var elemento = document.getElementById('target');
                  elemento.scrollIntoView(true); 
                
                  ")
            )
    
    # Obtener los usuarios de la base de datos
    sql_usuarios <- "SELECT * FROM usuarios"
    consulta_usuarios <- sqlInterpolate(conexion_base, sql_usuarios)
    usuarios_base<-dbGetQuery(conexion_base, consulta_usuarios)
    
    
    
    # Días para que la sesión expire
    cookie_caducidad<- 5
    
    # Función para guardas información de las sesiones en la base de datos
    add_sessionid_to_db <- function(user, sessionid, conn = conexion_base) {
        tibble(user = user, sessionid = sessionid, login_time = as.character(now())) %>%
            dbWriteTable(conn, "sessionids", ., append = TRUE)
    }
    
    # Función que regresa un data frame con la información del usuario y la sesión
    get_sessionids_from_db <- function(conn = conexion_base, expiry = cookie_caducidad) {
        dbReadTable(conn, "sessionids") %>%
            mutate(login_time = ymd_hms(login_time)) %>%
            as_tibble() %>%
            filter(login_time > now() - days(expiry))
    }

    #Modulo de la paquetería shinyauthr para el login de la página usando las credenciales creadas
    logout_init <- callModule(
        shinyauthr::logout,
        id = "logout",
        active = reactive(credenciales()$user_auth)
    )
    
    # Llama al modulo de login con las columnas de las contraseñas y usuarios de la base de datosy
    # crea las credenciales
    credenciales <- callModule(
        shinyauthr::login,
        id = "login",
        data = usuarios_base,
        user_col = user, 
        pwd_col = password,
        sessionid_col = sessionid, 
        cookie_getter = get_sessionids_from_db,
        cookie_setter = add_sessionid_to_db,
        sodium_hashed = TRUE,
        log_out = reactive(logout_init())
        )
    
    ############# Renderizado de la página en función del login y de las credenciales (usuarios, administradores)
    output$navbar<-
    renderUI(expr = if (credenciales()$user_auth==TRUE && user_data()$acceso==TRUE) {
        
        if(user_data()$permisos%in%"usuario"){
            #Aplicación si se registran como "usuario"
            navbarPage(title = "TopoApp", id="tabs", collapsible = TRUE,
                       tabPanel("Procesos",
                           h2(paste("Saludos ",user_data()$user)),    
                           tabsetPanel(
                               tabPanel("Ajuste Local",
                                        mainPanel(width = 12,
                                                  fluidRow(
                                                      column(width = 6, class="well",
                                                             h2("Sistema de origen"),
                                                             fileInput("archivoA", label = h3("Cargar datos"), accept = ".csv"),
                                                             uiOutput("tabla_origen")
                                                             
                                                      ),
                                                      column(width = 6, class="well",
                                                             h2("Sistema de destino"),
                                                             fileInput("archivoB", label = h3("Cargar datos"), accept = ".csv"),
                                                             uiOutput("tabla_destino")
                                                      )
                                                  ),
                                                  br(),
                                                  
                                                  fluidRow(
                                                      uiOutput("panel_emergente")
                                                  )
                                        )
                                        
                                        ),
                               tabPanel("UTM-Planas",
                                        fluidRow(
                                            column(10,
                                                   h2("Ajuste UTM-Planas"),
                                                   fileInput("archivoC", label = h3("Cargar datos"), accept = ".csv")
                                                   )
                                        ),
                                        fluidRow(
                                            column(12, class="well",
                                                   tabsetPanel(
                                                       tabPanel("Datos",
                                                                uiOutput("panel_utm_planas")
                                                       ),
                                                       tabPanel("Resultados",
                                                                uiOutput("panel_resultados")
                                                       )
                                                   )
                                            )
                                        )
                                        )
                           )
                       ),
                       tabPanel("Mapa",
                            sidebarPanel(width = 3,
                                         uiOutput("panel_mapa_local"),
                                         uiOutput("panel_mapa_UTM")
                                         ),
                            mainPanel(class="well",
                                       h2("Ubicación de puntos"),
                                       uiOutput("panel_mapa")
                            )  
                       )
                       
            )
            
        } else if(user_data()$permisos%in%"administrador"){
            # Aplicación si se ingresa como Administrador
            navbarPage(title = "TopoApp-Administrador", id="tabs1",
                       
                       tabPanel("Usuarios",
                                div(style="align-content: center;",
                                         column(width = 12,class="well",
                                                h3("Administrador de usuarios"),
                                                dataTableOutput("tabla_usuarios"),
                                                br(),
                                                br(),
                                                fluidRow(
                                                    column(width = 8, actionButton("acceso", label = "Acceso",  class = "btn-info")),
                                                    column(width = 2, actionButton("eliminar", label = "Eliminar",  class = "btn-info")),
                                                    column(width = 2, actionButton("nuevo", label = "Nuevo",  class = "btn-info"))
                                                )
                                                )
                                )
                                )
            )
            
        }
        
    })
    
    # Info de la sesión: tabla con los usuairios, sesiones, permisos, etc.
    user_data <- reactive({
        credenciales()$info
    })

    ## Conserva para ber la información del usuario
    output$user_table <- renderTable({
        # usar req para mostrar los resultados solo si: credentials()$user_auth es TRUE
        req(credenciales()$user_auth)
        user_data() %>%
            mutate(across(starts_with("login_time"), as.character))
    })
    
    # Renderiza la tabla de usuarios en la página de administrador
    usuarios<- reactiveValues()
    usuarios$datos<-usuarios_base
    
    output$tabla_usuarios<- renderDataTable({
        req(credenciales()$user_auth)
        
        usuarios_plot<- usuarios$datos[,c("user","email","permisos","nombre", "inicio","acceso")]
        Encoding(usuarios_plot$nombre)<-'UTF-8'
        datatable(usuarios_plot, options = list(
            language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json'),
            pageLength = 5,
            scrollX = TRUE,
            searching = FALSE
        ))
    })
    
    # Permíte a un nuevo usuario crearse una cuenta
    ##Pantalla emergenete para ingresar un Nuevo usuario
    observeEvent(input$registro,{
        showModal(
            modalDialog(
                title = "Únete",
                fluidRow(
                    column(width = 12,
                           textInput("usuarioNR", label = h4("Usuario")),
                           passwordInput("passwordNR",label = h4("Constraseña")),
                           textInput("emailNR", label = h4("Correo")),
                           textInput("nombreNR",label = h4("Nombre completo"))
                    )
                ),
                easyClose = FALSE,
                footer = tagList(
                    actionButton("cancelar","Cancelar"),
                    actionButton("guardarNR","Guardar")
                )
            )
        )
    })
    ##Ingresa un nuevo usuario a la base de datos sin antes aplicar algunas restricciones
    observeEvent(input$guardarNR,{
        if(nchar(input$usuarioNR)<1 || nchar(input$passwordNR)<1 ||
           nchar(input$emailNR)<1 || nchar(input$nombreNR)<1){
            showNotification(
                h4("Debes completar todos los campos"), 
                action = NULL, duration = 5, type = "warning")
        }else{
            if(length(grep("@",input$emailNR))==0){
                showNotification(
                    h4("Email debe tener un @"), 
                    action = NULL, duration = 5, type = "warning")
                
            } else {
                sql_usuarios <- "SELECT * FROM usuarios"
                usuarios_base<-dbGetQuery(conexion_base, sql_usuarios)
                
                if(length(grep(input$usuarioNR,usuarios_base$user))>0){
                    showNotification(
                        h4("El nombre de usuario ya existe, por favor prueba otro"), 
                        action = NULL, duration = 5, type = "warning")   
                }else{
                    
                    if(length(grep(input$emailNR,usuarios_base$email))>0){
                        showNotification(
                            h4("El correo ya existe, por favor ingresa con tu usuario registrado"), 
                            action = NULL, duration = 5, type = "warning")  
                    } else {
                        sql_nuevo_usuario <- 'INSERT INTO usuarios ("user","password",email,permisos,nombre,inicio, acceso) 
                    VALUES (?iduser,?idpassword,?idemail,?idpermiso,?idnombre,?idinicio,?idacceso);'
                        
                        consulta_usuarios <- sqlInterpolate(conexion_base, sql_nuevo_usuario, iduser= input$usuarioNR, 
                                                            idpassword=password_store(input$passwordNR),idemail=input$emailNR,
                                                            idpermiso="usuario",idnombre=input$nombreNR, idinicio=Sys.Date(),
                                                            idacceso="TRUE")
                        dbGetQuery(conexion_base, consulta_usuarios)
                        
                        removeModal()
                        refresh()
                    }
                    
                }
            }
        }
    })
    
    #Administrador CRUD
    ##Pantalla emergenete para ingresar un nuevo usuario
    observeEvent(input$nuevo,{
        showModal(
            modalDialog(
                title = "Nuevo usuario",
                fluidRow(
                    column(width = 12,
                           textInput("usuarioN", label = h4("Usuario")),
                           passwordInput("passwordN",label = h4("Constraseña")),
                           textInput("emailN", label = h4("Email")),
                           radioButtons("permisosN", label = h4("Permiso"),
                                        choiceNames = c("Usuario","Administrador"),
                                        choiceValues = c("usuario","administrador")
                                        ),
                           textInput("nombreN",label = h4("Nombres"))
                           )
                ),
                easyClose = FALSE,
                footer = tagList(
                    actionButton("cancelar","Cancelar"),
                    actionButton("guardarN","Guardar")
                )
            )
        )
    })
    ## Guarda los usuarios en la base de datos; antes aplica algunas restricciones para añdir un nuevo usuario
    observeEvent(input$guardarN,{
        if(nchar(input$usuarioN)<1 || nchar(input$passwordN)<1 ||
           nchar(input$emailN)<1 || nchar(input$nombreN)<1){
            showNotification(
                h4("Debes completar todos los campos"), 
                action = NULL, duration = 5, type = "warning")
        }else{
            if(length(grep("@",input$emailN))==0){
                showNotification(
                    h4("Email debe tener un @"), 
                    action = NULL, duration = 5, type = "warning")
                
            } else {
                sql_usuarios <- "SELECT * FROM usuarios"
                usuarios_base<-dbGetQuery(conexion_base, sql_usuarios)
                
                if(length(grep(input$usuarioN,usuarios_base$user))>0){
                    showNotification(
                        h4("El nombre de usuario ya existe"), 
                        action = NULL, duration = 5, type = "warning")   
                }else{
                    
                    if(length(grep(input$emailN,usuarios_base$email))>0){
                        showNotification(
                            h4("El correo ya existe"), 
                            action = NULL, duration = 5, type = "warning")  
                    } else {
                        sql_nuevo_usuario <- 'INSERT INTO usuarios ("user","password",email,permisos,nombre,inicio, acceso) 
                    VALUES (?iduser,?idpassword,?idemail,?idpermiso,?idnombre,?idinicio,?idacceso);'
                        
                        consulta_usuarios <- sqlInterpolate(conexion_base, sql_nuevo_usuario, iduser= input$usuarioN, 
                                                            idpassword=password_store(input$passwordN),idemail=input$emailN,
                                                            idpermiso=input$permisosN,idnombre=input$nombreN, idinicio=Sys.Date(),
                                                            idacceso="TRUE")
                        dbGetQuery(conexion_base, consulta_usuarios)
                        
                        removeModal()
                        
                        showNotification(
                            h4("Creación exitosa"), 
                            action = NULL, duration = 5, type = "message")
                        
                        sql_usuarios <- "SELECT * FROM usuarios"
                        usuarios$datos<-dbGetQuery(conexion_base, sql_usuarios)
                    }
                
                }
            }
        }
    })
    
    #Valores reactivos que guardaran info de las selecciones en la tabla
    selecciones_tabla<- reactiveValues()
    # Elimna un nuevo usuario
    observeEvent(input$eliminar,{
        if(length(input$tabla_usuarios_rows_selected)>0){
            showModal(
                modalDialog(title = "Borrar",
                            fluidPage(column(12,h3("Cuidado: Estás a punto de borrar usuarios de la base de datos"),style="color:red;")),
                            easyClose = FALSE,
                            size = "m",
                            footer = tagList(
                                actionButton("cancelar","Cancelar"),
                                actionButton("borrar_usuario","Eliminar")
                            ) 
                )
            )
        } else {
            showNotification(
                h4("Selecciona un renglón"), 
                action = NULL, duration = 5, type = "warning") 
        }
        
        
    })
    
    observeEvent(input$borrar_usuario,{
        selecciones_tabla$renglon<-input$tabla_usuarios_rows_selected
        usuarios_a_borrar<-usuarios$datos[selecciones_tabla$renglon,"user"]
        
        sql_borrar_usuario <- paste("DELETE FROM usuarios WHERE \"user\" IN (",
                                    paste0(sprintf("'%s'",usuarios_a_borrar),collapse = ","),")")
        
        dbGetQuery(conexion_base, sql_borrar_usuario)
        removeModal()
        
        showNotification(
            h4("Usuario eliminado con éxito"), 
            action = NULL, duration = 5, type = "message")
        
        sql_usuarios <- "SELECT * FROM usuarios"
        usuarios$datos<-dbGetQuery(conexion_base, sql_usuarios)
        
    })
    
    ## Cambia los accesos de usuarios
    observeEvent(input$acceso,{
        if(length(input$tabla_usuarios_rows_selected)>0){
            showModal(
                modalDialog(title = "Acceso",
                            fluidPage(
                                column(12,
                                       h3("Restringir de acceso a la aplicación"),
                                       radioButtons("acceso_usuario", label = h4("Acceso de usaurio"),
                                                    choiceNames = c("Restirngir acceso","Libre acceso"),
                                                    choiceValues = c(FALSE, TRUE)
                                                    )
                                       )
                                ),
                            easyClose = FALSE,
                            size = "m",
                            footer = tagList(
                                actionButton("cancelar","Cancelar"),
                                actionButton("cambiar_acceso","Cambiar")
                            ) 
                )
            )
        } else {
            showNotification(
                h4("Selecciona un renglón"), 
                action = NULL, duration = 5, type = "warning") 
        }  
    })
    
    observeEvent(input$cambiar_acceso,{
        selecciones_tabla$renglon<-input$tabla_usuarios_rows_selected
        usuarios_con_acceso<-usuarios$datos[selecciones_tabla$renglon,"user"]
        
        sql_acceso <- paste0("UPDATE usuarios SET \"acceso\" ='",input$acceso_usuario,"'  WHERE \"user\" IN (",
                                    paste0(sprintf("'%s'",usuarios_con_acceso),collapse = ","),")")
        
        dbGetQuery(conexion_base, sql_acceso)
        removeModal()
        
        showNotification(
            h4("Cambio exitoso"), 
            action = NULL, duration = 5, type = "message")
        
        sql_usuarios <- "SELECT * FROM usuarios"
        usuarios$datos<-dbGetQuery(conexion_base, sql_usuarios)
        
    })
    
    #boton general que cierra las ventanas emergentes si se selecciona el botón cancelar
    observeEvent(input$cancelar,{
        removeModal()
    })
    
    
################################Procesos ###########################
###############proceso local#######################    
    # Genera las tablas cuando estas son cargadas a la aplicación
    ## Valores reactivos
    datos<- reactiveValues()
    ##tabla izquirda
    output$tabla_origen<-
        renderUI(expr = if (!is.null(input$archivoA)) {
            dataTableOutput("tabla_puntos_origen", )
        } else {
            NULL
        })
    output$tabla_puntos_origen<-renderDataTable({
        req(credenciales()$user_auth)
        datos$puntos_origen<- read.csv(input$archivoA$datapath, sep = ",", header = FALSE)
        datatable(datos$puntos_origen, options = list(
            language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json'),
            pageLength = 5,
            scrollX = TRUE,
            searching = FALSE
        ))
    })
    ##otra tabla
    output$tabla_destino<-
        renderUI(expr = if (!is.null(input$archivoB)) {
            dataTableOutput("tabla_puntos_destino")
        } else {
            NULL
        })
    
    output$tabla_puntos_destino<-renderDataTable({
        req(credenciales()$user_auth)
        datos$puntos_destino<- read.csv(input$archivoB$datapath, sep = ",", header = FALSE)
        datatable(datos$puntos_destino, options = list(
            language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json'),
            pageLength = 5,
            scrollX = TRUE,
            searching = FALSE
        ))
    })
    #Apartado cálculos, resultados y modelo
    output$panel_emergente<-
        renderUI(expr = if (!is.null(input$archivoB) && !is.null(input$archivoB)) {
            fluidPage(
                sidebarPanel(width = 3,
                             h2("Ajustes"),
                             checkboxGroupInput("proceso_local", label = h3("Proceso local"),
                                                choiceNames = c("Horizontal","Vertical"),
                                                choiceValues = c("horizontal","vertical"),
                                                selected = "horizontal"
                             ),
                             uiOutput("botones_coordenadas")
                ),
                mainPanel(class="well", id="target",
                    tabsetPanel(type = "tabs",
                                tabPanel("Puntos de control",
                                         fluidRow(
                                             h4("Selecciona el nombre de las columnas que hacen match"),
                                             column(3, selectInput("col_name_origen", label = h4("Columna de nombres"),
                                                                   choices = names(datos$puntos_origen))),
                                             column(3, selectInput("col_name_destino", label = h4("Columna de nombres"),
                                                                   choices = names(datos$puntos_destino))),
                                             column(3,
                                                    h4("Crear puntos de control"),
                                                    actionButton("puntos_control","Crear", class = "btn-danger")
                                             ),
                                             dataTableOutput("tabla_puntos_match")
                                             
                                         )
                                ),
                                tabPanel("Modelo", verbatimTextOutput("modelo"),
                                         downloadLink('descarga_modelo', 'Descargar datos del modelo')
                                         ),
                                tabPanel("Nuevos cálculos",
                                         fluidRow(
                                             column(9,
                                                    fileInput("archivo_corregir", label = h3("Cargar datos para correción"), accept = ".csv"),
                                                    dataTableOutput("datos_a_corregir")
                                                    ),
                                             column(3, uiOutput("col_datos_correccion"))
                                         ),
                                         br(),
                                         br(),
                                         uiOutput("panel_emergente_resultados_local")
                                         )
                    )
                )
            )
            
        } else {
            NULL
        })
    ## Tabla para matches en función del nombre que escoga el usuario
    
    observeEvent(input$puntos_control,{
        
        #ID  para quue al momento de hacer el match no se ordenen
        origen<-datos$puntos_origen
        origen$id  <- 1:nrow(origen)
        matchP<-merge(origen, datos$puntos_destino, by.x=input$col_name_origen, by.y=input$col_name_destino)
        matchP<-matchP[order(matchP$id),]
        datos$datos_match<-matchP%>%select(-id)
        
        output$tabla_puntos_match<-renderDataTable({
            req(credenciales()$user_auth)  
            datatable(datos$datos_match, options = list(
                language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json'),
                pageLength = 5,
                scrollX = TRUE,
                searching = FALSE
            ))
        })
        
    })
    
    ##Panel para pedir coordendas para construir la matriz
    output$botones_coordenadas<-renderUI(
        expr = if (is.null(datos$datos_match) || is.null(input$proceso_local) || nrow(datos$datos_match)<1) {
            NULL
        } else {
            if(length(input$proceso_local)>1){
                
                div(h4("Selecciona las columnas para construir la matriz"), 
                         selectInput("col_N", label = "Columna N(Norte)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[2]),
                         selectInput("col_E", label = "Columna E(Este)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[3]),
                         selectInput("col_h", label = "Columna h(Altura elipsoidal)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[4]),
                         selectInput("col_Y", label ="Columna Y(Norte)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[5]),
                         selectInput("col_X", label = "Columna X(Este)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[6]),
                         selectInput("col_H", label = "Columna H(Altura ortométrica)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[7]),
                         actionButton("inicio_ajuste_local",label = "Iniciar", class = "btn-danger")
                )
            }else if(input$proceso_local=="horizontal"){
                div(h4("Selecciona las columnas para construir la matriz"), 
                         selectInput("col_N", label = "Columna N(Norte)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[2]),
                         selectInput("col_E", label = "Columna E(Este)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[3]),
                         selectInput("col_Y", label = "Columna Y(Norte)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[4]),
                         selectInput("col_X", label = "Columna X(Este)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[5]),
                         actionButton("inicio_ajuste_local",label = "Iniciar", class = "btn-danger")
                )
                
            } else if (input$proceso_local=="vertical"){
                div(h4("Selecciona las columnas para construir la matriz"),
                         selectInput("col_N", label = "Columna N(Norte)",
                                     choices = names(datos$datos_match), selected=names(datos$datos_match)[2]),
                         selectInput("col_E", label = "Columna E(Este)",
                                     choices = names(datos$datos_match), selected=names(datos$datos_match)[3]),
                         selectInput("col_h", label = "Columna h(Altura elipsoidal)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[4]),
                         selectInput("col_H", label = "Columna H(Altura ortométrica)",
                                               choices = names(datos$datos_match), selected=names(datos$datos_match)[5]),
                         actionButton("inicio_ajuste_local",label = "Iniciar",class = "btn-danger"),
                )
            } else {
                NULL
            }
        }
        )
    
    ######## Cálculo Ajuste local (Obtener modelo)################
    
    observeEvent(input$aceptar_proceso_local,{
        if(length(input$proceso_local)>1){
            
            matriz_coordenadas_H<-datos$datos_match[,c(input$col_X, input$col_Y,input$col_E, input$col_N)]
            matriz_coordenadas_V<-datos$datos_match[,c(input$col_N, input$col_E,input$col_H, input$col_h)]
            
            # Restricciones
            if(class(matriz_coordenadas_H[,input$col_X])!="numeric" || class(matriz_coordenadas_H[,input$col_Y])!="numeric" ||
               class(matriz_coordenadas_H[,input$col_E])!="numeric" || class(matriz_coordenadas_H[,input$col_N])!="numeric" ||
               class(matriz_coordenadas_V[,input$col_N])!="numeric" || class(matriz_coordenadas_V[,input$col_E])!="numeric" ||
               class(matriz_coordenadas_V[,input$col_H])!="numeric" || class(matriz_coordenadas_V[,input$col_h])!="numeric"
               ){
                
                showModal(modalDialog(title = "Error",
                                      h4("Error: Tus columnas no son de tipo número: a)Asugurate de que tus columnas sean de tipo número
                             b) Asegurate de que tus tablas no tengan un encabezado"
                                      ),
                                      size = "m",
                                      easyClose = TRUE
                )
                )
                return()
            }
            
            matriz_resultados_H<-crear_X(matriz_coordenadas_H[,c(3,4,1,2)])
            names(matriz_resultados_H)<-c("Matriz_A","Parámetros","Vector_L","Vector_residuos")
            
            matriz_resultados_V<-crear_B(matriz_coordenadas_V)
            names(matriz_resultados_V)<-c("Parámetros","Vector_Y","Vector_B")
            
            matriz_resultados_HV<-list(matriz_resultados_H, matriz_resultados_V)
            names(matriz_resultados_HV)<-c("Horizontal","Vertical")
            
            datos$matrices<-matriz_resultados_HV
            
        } else if (input$proceso_local=="horizontal"){
            matriz_coordenadas<-datos$datos_match[,c(input$col_X, input$col_Y,input$col_E, input$col_N)]
            
            # Restricciones
            if(class(matriz_coordenadas[,input$col_X])!="numeric" || class(matriz_coordenadas[,input$col_Y])!="numeric" ||
               class(matriz_coordenadas[,input$col_E])!="numeric" || class(matriz_coordenadas[,input$col_N])!="numeric"){
                
                showModal(modalDialog(title = "Error",
                                      h4("Error: Tus columnas no son de tipo número: a)Asugurate de que tus columnas sean de tipo número
                             b) Asegurate de que tus tablas no tengan un encabezado"
                                      ),
                                      size = "m",
                                      easyClose = TRUE
                )
                )
                return()
            }
            
            matriz_resultados<-crear_X(matriz_coordenadas[,c(3,4,1,2)])
            names(matriz_resultados)<-c("Matriz_A","Parámetros","Vector_L","Vector_residuos") ##A,X, L, V
            datos$matrices<-matriz_resultados
        } else if(input$proceso_local=="vertical"){
            matriz_coordenadas<-datos$datos_match[,c(input$col_N, input$col_E,input$col_H, input$col_h)]
            
            # Restricciones
            if(class(matriz_coordenadas[,input$col_N])!="numeric" || class(matriz_coordenadas[,input$col_E])!="numeric" ||
               class(matriz_coordenadas[,input$col_H])!="numeric" || class(matriz_coordenadas[,input$col_h])!="numeric"){
                
                showModal(modalDialog(title = "Error",
                          h4("Error: Tus columnas no son de tipo número: a)Asugurate de que tus columnas sean de tipo número
                             b) Asegurate de que tus tablas no tengan un encabezado"
                             ),
                          size = "m",
                          easyClose = TRUE
                          )
                          )
                return()
            }
            
            matriz_resultados<-crear_B(matriz_coordenadas)
            names(matriz_resultados)<-c("Parámetros","Vector_Y","Vector_B")
            datos$matrices<-matriz_resultados
        }else{
            NULL
        }
        
        removeModal()
    })
    
    output$modelo<-renderPrint({datos$matrices})
    
    ## Descarga los resultados del modelo
    output$descarga_modelo <- downloadHandler(
        
        filename = function() {"modelo.txt"},
        content = function(file) {
            
            sink(file,append = TRUE)
            print(datos$matrices)
            sink()
        })
    
    ######## Cálculo Ajuste local (Correción de valores)################
    # Tabla que se crea al momento de cargar los datos a corregir
    output$datos_a_corregir<-renderDataTable({
        if(is.null(input$archivo_corregir)){
            NULL
        } else {
            datos$datos_corregir<-read.csv(input$archivo_corregir$datapath, sep = ",", header = FALSE)
            datatable(datos$datos_corregir, options = list(
                language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json'),
                pageLength = 5,
                scrollX = TRUE,
                searching = FALSE
            ))
        }
    })
    # Panel emergete para seleccionar las columnas para hacer el cálculo de correcciones
    output$col_datos_correccion<-renderUI(
        expr = if (is.null(input$archivo_corregir) || length(input$proceso_local)==0 || is.null(datos$matrices)) {
            NULL
        } else {
            if(length(input$proceso_local)>1){
                
                fluidRow(h4("Selecciona las columnas para corrección"), class="well",
                         radioButtons("panel_correccion_multiple", label = "Selecciona la correción",
                                      choiceNames = c("Horizontal","Vertical"),
                                      choiceValues = c("horizontal","vertical")
                                      ),
                         uiOutput("panel_correccion_nuevo")
                )
            }else if(input$proceso_local=="horizontal"){
                fluidRow(h4("Selecciona las columnas para corrección"),
                         selectInput("col_nom_local", label = "Columna Punto",
                                     choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[1]),
                         selectInput("col_Y_c", label = "Columna Y(Norte)",
                                     choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[2]),
                         selectInput("col_X_c", label = "Columna X(Este)",
                                     choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[3]),
                         actionButton("inicio_correccion_local",label = "Iniciar", class = "btn-danger")
                )
                
            } else if (input$proceso_local=="vertical"){
                fluidRow(h4("Selecciona las columnas para corrección"),
                         selectInput("col_nom_local", label = "Columna Punto",
                                     choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[1]),
                         selectInput("col_N_c", label = "Columna N(Norte)",
                                     choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[2]),
                         selectInput("col_E_c", label = "Columna E(Este)",
                                     choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[3]),
                         selectInput("col_h_c", label = "Columna h(Altura elipsoidal)",
                                     choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[4]),
                         actionButton("inicio_correccion_local",label = "Iniciar", class = "btn-danger")
                )
            } else {
                NULL
            }
        }
    )
    
    ### Panel emergente para corregir datos si el modelo fue múltiple
    output$panel_correccion_nuevo<- renderUI(
        expr = if(input$panel_correccion_multiple=="vertical"){
            fluidRow(
                column(10,
                       selectInput("col_nom_local", label = "Columna Punto",
                                   choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[1]),
                       selectInput("col_N_c", label = "Columna N(Norte)",
                                   choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[2]),
                       selectInput("col_E_c", label = "Columna E(Este)",
                                   choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[3]),
                       selectInput("col_h_c", label = "Columna h(Altura elipsoidal)",
                                   choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[4]),
                       actionButton("inicio_correccion_local",label = "Iniciar", class = "btn-danger")
                       )
            )
        } else if(input$panel_correccion_multiple=="horizontal"){
            fluidRow(
                column(10,
                       selectInput("col_nom_local", label = "Columna Punto",
                                   choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[1]),
                       selectInput("col_Y_c", label = "Columna Y(Norte)",
                                   choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[2]),
                       selectInput("col_X_c", label = "Columna X(Este)",
                                   choices = names(datos$datos_corregir), selected=names(datos$datos_corregir)[3]),
                       actionButton("inicio_correccion_local",label = "Iniciar", class = "btn-danger")
                       )
            )
        }
    )
    
    ######## Cálculo Ajuste local (Obtener nuevos valores corregidos)################
    
    observeEvent(input$aceptar_calculo_local,{
        if(length(input$proceso_local)>1){
            if(input$panel_correccion_multiple=="horizontal"){
                ## Restricción para evitar que los calculos se hagan con matrices no correspondientes
                if(length(datos$matrices)==4){
                    showNotification(
                        h4("Error: Cambia el tipo de proceso o haz un nuevo cálculo"), 
                        action = NULL, duration = 5, type = "warning")
                    return()
                }
                matriz_corregir<-datos$datos_corregir[,c(input$col_X_c, input$col_Y_c)]
                matriz_resultados_correccion<-resultados(matriz_corregir,datos$matrices[[1]][[2]])
                datos$datos_corregidos<-matriz_resultados_correccion[,c(2,1)]
                
                removeModal()
                
            } else {
                if(length(datos$matrices)==4){
                    showNotification(
                        h4("Error: Cambia el tipo de proceso o haz un nuevo cálculo"), 
                        action = NULL, duration = 5, type = "warning")
                    return()
                }
                matriz_corregir<-datos$datos_corregir[,c( input$col_N_c, input$col_E_c,input$col_h_c)]
                matriz_resultados_correccion<-resultados_H(matriz_corregir,datos$matrices[[2]][[3]])
                datos$datos_corregidos<-matriz_resultados_correccion[[2]] 
                
                removeModal()
            }
            
        } else if (input$proceso_local=="horizontal"){
            ## Restricción para evitar que los calculos se hagan con matrices no correspondientes
            if(length(datos$matrices)==2 || length(which("V" %in% names(datos$matrices)))==0){
                showNotification(
                    h4("Error: Cambia el tipo de proceso o haz un nuevo cálculo"), 
                    action = NULL, duration = 5, type = "warning")
                return()
            }
            matriz_corregir<-datos$datos_corregir[,c(input$col_X_c, input$col_Y_c)]
            matriz_resultados_correccion<-resultados(matriz_corregir,datos$matrices[[2]])
            datos$datos_corregidos<-matriz_resultados_correccion[,c(2,1)]
            
            removeModal()
        } else if(input$proceso_local=="vertical"){
            ## Restricción para evitar que los calculos se hagan con matrices no correspondientes
            if(length(datos$matrices)==2 || length(which("B" %in% names(datos$matrices)))==0){
                showNotification(
                    h4("Error: Cambia el tipo de proceso o haz un nuevo cálculo"), 
                    action = NULL, duration = 5, type = "warning")
                return()
            }
            matriz_corregir<-datos$datos_corregir[,c( input$col_N_c, input$col_E_c,input$col_h_c)]
            matriz_resultados_correccion<-resultados_H(matriz_corregir,datos$matrices[[3]])
            datos$datos_corregidos<-matriz_resultados_correccion[[2]]
            
            removeModal()
        }else{
            NULL
        }
    })
    
    
    output$texto_corregido<-renderPrint({cbind(datos$datos_corregir[,input$col_nom_local],datos$datos_corregidos)})
    
    ## Genera el panel inferior al momento de calcular nuevos datos
    output$panel_emergente_resultados_local<- renderUI(
        expr = if(is.null(datos$datos_corregidos)){
            NULL
        } else {
            fluidRow(
                column(10, verbatimTextOutput("texto_corregido"),
                       downloadLink('descarga_local_resultado', 'Descargar resultados')
                )
                
            ) 
        }
    )
    
    ## Descarga los nuevos datos del proceso local
    output$descarga_local_resultado<-downloadHandler(
        filename="resultados_proceso_local.csv",
        content=function(file){
            csv_resultados_local<- cbind(datos$datos_corregir[,input$col_nom_local],datos$datos_corregidos)
            
            write.csv(csv_resultados_local,file, row.names = FALSE)
        }
    )
    
############proceso UTM-PLANAS######################
    ## Genera el panel de procesos al momento de cargar los datos
    output$panel_utm_planas<- renderUI(expr = if (!is.null(input$archivoC)) {
            
        fluidPage(
            sidebarPanel(width = 3,
                         selectInput("crs", "CRS",
                                     list("WGS84-UTM-12 Norte" = 32612, "WGS84-UTM-13 Norte" = 32613, "WGS84-UTM-14 Norte" = 32614,
                                          "WGS84-UTM-15 Norte" = 32615, "WGS84-UTM-16 Norte" = 32616, "WGS84-UTM-17 Sur" = 32717,
                                          "WGS84-UTM-18 Sur" = 32718, "WGS84-UTM-19 Sur" = 32719, "WGS84-UTM-20 Sur" = 32720
                                     )
                         ),
                         selectInput("col_nombre_UTM", label = "Columna Punto",
                                     choices = names(datos$puntos_coordenadasUTM), selected=names(datos$puntos_coordenadasUTM)[1]),
                         selectInput("col_N_UTM", label = "Columna N(Norte)",
                                     choices = names(datos$puntos_coordenadasUTM), selected=names(datos$puntos_coordenadasUTM)[2]),
                         selectInput("col_E_UTM", label = "Columna E(Este)",
                                     choices = names(datos$puntos_coordenadasUTM), selected=names(datos$puntos_coordenadasUTM)[3]),
                         selectInput("col_h_UTM", label = "Columna h(Altura elipsoidal)",
                                     choices = names(datos$puntos_coordenadasUTM), selected=names(datos$puntos_coordenadasUTM)[4]),
                         actionButton("inicio_correccion_UTM",label = "Iniciar", class = "btn-info")
            ),
            mainPanel( class="well",
                       h4("Selecciona un renglón en la tabla (Punto pivote)"),
                       dataTableOutput("tabla_inicio_utm")
            )
        )   
    } else {
        NULL
    })
    
    ## Geera la tabla al momento de cargar los datos
    output$tabla_inicio_utm<- renderDataTable({
        req(credenciales()$user_auth)
        datos$puntos_coordenadasUTM<-read.csv(input$archivoC$datapath, sep = ",", header = FALSE)
        datatable(datos$puntos_coordenadasUTM, options = list(
            language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json'),
            pageLength = 5,
            scrollX = TRUE,
            searching = FALSE
        ), selection='single')
    })
    
    ####### Realiza el ajuste de UTM planas
    observeEvent(input$aceptar_proceso_utm,{
        
        
        # Restricciones
        if(length(input$tabla_inicio_utm_rows_selected)==0){
            showNotification(
                h4("Selecciona un renglón correspondiente con el punto pivote"), 
                action = NULL, duration = 5, type = "warning")
            
            removeModal()
            return()
        }
        
        datos$datos_utm_ordenados<-rbind(datos$puntos_coordenadasUTM[input$tabla_inicio_utm_rows_selected,], datos$puntos_coordenadasUTM[-input$tabla_inicio_utm_rows_selected,])
        datos_utm<-datos$datos_utm_ordenados[,c(input$col_nombre_UTM,input$col_E_UTM, input$col_N_UTM,input$col_h_UTM)]
        names(datos_utm)<-c("Punto","E","N","q")
        
        #Restricciones
        if(class(datos_utm$N)!="numeric" && class(datos_utm$E)!="numeric" &&
           class(datos_utm$q)!="numeric"){
            
            showModal(modalDialog(title = "Error",
                                  h4("Error: Tus columnas no son de tipo número: a)Asugurate de que tus columnas sean de tipo número
                             b) Asegurate de que tus tablas no tengan un encabezado"
                                  ),
                                  size = "m",
                                  easyClose = TRUE
            ))
            return()
        }
        
        
        datos$correccion_utm<-funcion_UTM_planas(datos_utm[,-1], as.numeric(input$crs))%>%select(-Acimut_rad,-acimut_grad)
        
        removeModal()
    })
    
    output$tabla_utm_corregido<- renderDataTable({
        
        datos$datos_utm_ordenados<-rbind(datos$puntos_coordenadasUTM[input$tabla_inicio_utm_rows_selected,], datos$puntos_coordenadasUTM[-input$tabla_inicio_utm_rows_selected,])
        datos_utm<-datos$datos_utm_ordenados[,c(input$col_nombre_UTM,input$col_E_UTM, input$col_N_UTM,input$col_h_UTM)]
        names(datos_utm)<-c("Punto","N","E","h")
        datos_tabla_UTM<-cbind(datos_utm,datos$correccion_utm)
        
        datatable(datos_tabla_UTM, options = list(
            language = list(url = '//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json'),
            pageLength = 5,
            scrollX = TRUE,
            searching = FALSE
        ))
    })
    ## Genera un nuevo panel donde se muestran los resultaos del ajuste
    output$panel_resultados<-renderUI(expr = if(!is.null(datos$correccion_utm)){
        fluidRow(
            column( 10, class="well",
                h3("Resultados"),
                dataTableOutput("tabla_utm_corregido"),
                downloadLink('descarga_utm_correccion', 'Descargar resultados')
            )
        )
    } else{
        NULL
    })
    
    ## Descarga los resultados del proceso UTM-planas
    output$descarga_utm_correccion<-downloadHandler(
        filename="resultados_proceso_UTM-planas.csv",
        content=function(file){
            csv_resultados_local<- cbind(datos$datos_utm_ordenados,datos$correccion_utm[,c("k_esc","K_ele", "K_com", "N_top","E_top")])
            
            write.csv(csv_resultados_local,file, row.names = FALSE, col.names=FALSE)
        }
    )
    
    
    ##################### Mapa ########################
    ## Genera la ventana del mapa cuando los datos fueron creados correctamente
    output$panel_mapa<-renderUI(
        expr = if(is.null(datos$mapa_datos_input) && is.null(datos$mapa_datos_resultados)){
            NULL
        } else {
            leafletOutput("mapa")
        }
    )
    ## Crea el mapa
    output$mapa<- renderLeaflet({
        mapa<-leaflet() %>%
            addProviderTiles(providers$OpenStreetMap.Mapnik, group = "OpenStreetMap.Mapnik") %>%
            addProviderTiles(providers$Esri.WorldImagery, group = "Esri.WorldImagery") %>%
            addLayersControl(
                baseGroups = c("OpenStreetMap.Mapnik","Esri.WorldImagery"),
                overlayGroups = c("Puntos de control", "Resultados")
            )%>%
            addCircleMarkers(data=datos$mapa_datos_input, color = "red", group = "Puntos de control",
                             label = ~as.character(datos$mapa_datos_input[,3]))%>%
            addCircleMarkers(data=datos$mapa_datos_resultados, color="blue", group = "Resultados",
                             label = ~as.character(datos$mapa_datos_input[,3]))
        
        mapa
    })
    
    # Genera el panel para crear el mapa de los datos de la corrección UTM
    output$panel_mapa_UTM<- renderUI(
        expr = if(is.null(datos$correccion_utm)){
            NULL
        } else {
            fluidRow(
                column(10,
                       h3("Datos del proceso UTM-Planas"),
                       selectInput("crs_mapa_UTM", "Selecciona un sistema de referencia; Asegurate de que se al mismo que usaste en los procesos",
                                   list("WGS84-UTM-12 Norte" = 32612, "WGS84-UTM-13 Norte" = 32613, "WGS84-UTM-14 Norte" = 32614,
                                        "WGS84-UTM-15 Norte" = 32615, "WGS84-UTM-16 Norte" = 32616, "WGS84-UTM-17 Sur" = 32717,
                                        "WGS84-UTM-18 Sur" = 32718, "WGS84-UTM-19 Sur" = 32719, "WGS84-UTM-20 Sur" = 32720
                                   )
                       ),
                       actionButton("crear_mapa_utm","Iniciar")
                )
            )
        }
    )
    
    observeEvent(input$crear_mapa_utm,{
        datos_utm<-datos$datos_utm_ordenados[,c(input$col_nombre_UTM,input$col_E_UTM, input$col_N_UTM,input$col_h_UTM)]
        ## Restricciones
        if(class(datos_utm[,2])!="numeric" || class(datos_utm[,3])!="numeric" || class(datos_utm[,4])!="numeric"){
            showNotification(
                h4("Las columnas que seleccionaste en el proceso no son numéricas"), 
                action = NULL, duration = 5, type = "warning")
            return()
        }
        datos_mapa_UTM<-cbind(datos_utm,datos$correccion_utm)
        sf_datos_locales<-st_as_sf(datos_mapa_UTM, coords=c("Long","Lat"), crs=4326)
        
        datos$mapa_datos_input<-sf_datos_locales
        
        sf_datos_resultados<-st_as_sf(datos_mapa_UTM, coords=c("E_top","N_top"), crs=as.numeric(input$crs_mapa_UTM))
        sf_datos_resultados<-st_transform(sf_datos_resultados, 4326)
        
        datos$mapa_datos_resultados<-sf_datos_resultados
        
    })
    
    # Genera el panel para crear el mapa de los datos de la corrección local
    output$panel_mapa_local<- renderUI(
        expr = if(is.null(datos$datos_corregidos)){
            NULL
        } else {
            fluidRow(
                column(10,
                       h3("Datos del proceso local"),
                       selectInput("crs_mapa_local", "Selecciona un sistema de referencia de coordenadas",
                                   list("WGS84-UTM-12 Norte" = 32612, "WGS84-UTM-13 Norte" = 32613, "WGS84-UTM-14 Norte" = 32614,
                                        "WGS84-UTM-15 Norte" = 32615, "WGS84-UTM-16 Norte" = 32616, "WGS84-UTM-17 Sur" = 32717,
                                        "WGS84-UTM-18 Sur" = 32718, "WGS84-UTM-19 Sur" = 32719, "WGS84-UTM-20 Sur" = 32720
                                   )
                       ),
                       actionButton("crear_mapa_local","Iniciar")
                       )
            )
        }
    )
    
    observeEvent(input$crear_mapa_local,{
        resultados_utm<-cbind(datos$datos_corregir[,input$col_nom_local],datos$datos_corregidos)
        resultados_utm[,2]<-as.numeric(resultados_utm[,2])
        resultados_utm[,3]<-as.numeric(resultados_utm[,3])
        resultados_utm_sf<-resultados_utm %>%as.data.frame%>%st_as_sf(coords=c(2,3), crs= as.numeric(input$crs_mapa_local))
        resultados_utm_sf<- st_transform(resultados_utm_sf, 4326)
        
        datos$mapa_datos_resultados<-resultados_utm_sf
        
        ## Restricciones
        if(class(datos$datos_match[,c(input$col_E)])!="numeric" || class(datos$datos_match[,c(input$col_N)])!="numeric" ){
            showNotification(
                h4("Las columnas que seleccionaste en el proceso no son numéricas"), 
                action = NULL, duration = 5, type = "warning")
            return()
        }
        match_sf<-st_as_sf(datos$datos_match, coords=c(input$col_E, input$col_N), crs=as.numeric(input$crs_mapa_local))
        match_sf<-st_transform(match_sf, 4326)
        
        datos$mapa_datos_input<-match_sf
    })
    
    
    
    ###########################################Guarda los datos en la base de datos####################
    
    observeEvent(input$inicio_ajuste_local,{
        showModal(modalDialog(title="Guardar",
                              fluidRow(h4("Asgurate de haber elegido los nombres de las coordenadas correctamente antes de continuar")),
                              footer = tagList(
                                  actionButton("cancelar","Cancelar"),
                                  actionButton("aceptar_proceso_local","Continuar")
                              )
                              ))
    })
    
    observeEvent(input$aceptar_proceso_local,{
        if(length(input$proceso_local)>1){
            guardar_match<-datos$datos_match[,c(input$col_N, input$col_E, input$col_h,input$col_Y, input$col_X,input$col_H)]
            names(guardar_match)<-c("n","e","he","y","x","ho")
            ### Restricciones
            if(class(guardar_match$n)!="numeric" || class(guardar_match$e)!="numeric" || class(guardar_match$he)!="numeric" ||
               class(guardar_match$y)!="numeric" || class(guardar_match$x)!="numeric" || class(guardar_match$ho)!="numeric"){
                showNotification(
                    h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                    action = NULL, duration = 5, type = "warning")
                return()
            }
        }else if(input$proceso_local%in%"horizontal"){
            guardar_match<-datos$datos_match[,c(input$col_N, input$col_E, input$col_Y, input$col_X)]
            names(guardar_match)<-c("n","e","y","x")
            ### Restricciones
            if(class(guardar_match$n)!="numeric" || class(guardar_match$e)!="numeric" || 
               class(guardar_match$y)!="numeric" || class(guardar_match$x)!="numeric"){
                showNotification(
                    h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                    action = NULL, duration = 5, type = "warning")
                return()
            }
        } else if(input$proceso_local%in%"vertical"){
            guardar_match<-datos$datos_match[,c(input$col_N, input$col_E, input$col_h, input$col_H)]
            names(guardar_match)<-c("n","e","he","ho")
            ### Restricciones
            if(class(guardar_match$n)!="numeric" || class(guardar_match$e)!="numeric" || class(guardar_match$he)!="numeric" ||
               class(guardar_match$ho)!="numeric"){
                showNotification(
                    h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                    action = NULL, duration = 5, type = "warning")
                return()
            }
        } else {
            NULL
        }
        
        user_id<-data.frame(user_id=rep(user_data()$user_id,length(datos$datos_match[,1])))
        fecha<-data.frame(fecha=rep(Sys.Date(),length(datos$datos_match[,1])))
        ##Restricciones
        if(class(datos$datos_match[,1])!="character"){
            showNotification(
                h4("Las columnas nombre de puto no es caracter: No se guardaran los datos"), 
                action = NULL, duration = 5, type = "warning")
            return()
        }
        punto<-datos$datos_match[,1]
        
        guardar_datos<-cbind(user_id,fecha,punto,guardar_match)
        guardar_datos$fecha<-as.character(guardar_datos$fecha)
        query_agregar<-sqlAppendTable(conexion_base, 'puntos_control_local',guardar_datos )
        
        dbGetQuery(conexion_base,query_agregar)
        
        
        showNotification(
            h4("Datos guardados con éxito"), 
            action = NULL, duration = 5, type = "message")
        
    })
    
    ## Guarda la tabla de carga de nuevos datos
    observeEvent(input$inicio_correccion_local,{
        showModal(modalDialog(title="Guardar",
                              fluidRow(h4("Asgurate de haber elegido los nombres de las coordenadas correctamente antes de continuar")),
                              footer = tagList(
                                  actionButton("cancelar","Cancelar"),
                                  actionButton("aceptar_calculo_local","Continuar")
                              )
        ))
    })
    
    observeEvent(input$aceptar_calculo_local,{
        if(length(input$proceso_local)>1){
            
            if(input$panel_correccion_multiple%in%"horizontal"){
                guardar_match<-datos$datos_corregir[,c(input$col_Y_c, input$col_X_c)]
                names(guardar_match)<-c("y","x")
                ### Restricciones
                if(class(guardar_match$y)!="numeric" || class(guardar_match$x)!="numeric"){
                    showNotification(
                        h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                        action = NULL, duration = 5, type = "warning")
                    return()
                }
                
            }else if(input$panel_correccion_multiple%in%"vertical"){
                guardar_match<-datos$datos_corregir[,c(input$col_N_c, input$col_E_c, input$col_h_c)]
                names(guardar_match)<-c("n","e","he") 
                ### Restricciones
                if(class(guardar_match$n)!="numeric" || class(guardar_match$e)!="numeric" || 
                   class(guardar_match$he)!="numeric" ){
                    showNotification(
                        h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                        action = NULL, duration = 5, type = "warning")
                    return()
                }
            } else { 
                NULL
                }
        }else if(input$proceso_local%in%"horizontal"){
            guardar_match<-datos$datos_corregir[,c(input$col_Y_c, input$col_X_c)]
            names(guardar_match)<-c("y","x")
            ### Restricciones
            if(class(guardar_match$y)!="numeric" || class(guardar_match$x)!="numeric"){
                showNotification(
                    h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                    action = NULL, duration = 5, type = "warning")
                return()
            }
        } else if(input$proceso_local%in%"vertical"){
            guardar_match<-datos$datos_corregir[,c(input$col_N_c, input$col_E_c, input$col_h_c)]
            names(guardar_match)<-c("n","e","he")
            ### Restricciones
            if(class(guardar_match$n)!="numeric" || class(guardar_match$e)!="numeric" || 
               class(guardar_match$he)!="numeric"){
                showNotification(
                    h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                    action = NULL, duration = 5, type = "warning")
                return()
            }
        } else {
            NULL
        }
        
        user_id<-data.frame(user_id=rep(user_data()$user_id,length(datos$datos_corregir[,1])))
        fecha<-data.frame(fecha=rep(Sys.Date(),length(datos$datos_corregir[,1])))
        ##Restricciones
        if(class(datos$datos_corregir[,c(input$col_nom_local)])!="character"){
            showNotification(
                h4("La columna nombre de puto no es caracter: No se guardaran los datos"), 
                action = NULL, duration = 5, type = "warning")
            return()
        }
        punto<-datos$datos_corregir[,c(input$col_nom_local)]
        
        guardar_datos<-cbind(user_id,fecha,punto,guardar_match)
        guardar_datos$fecha<-as.character(guardar_datos$fecha)
        query_agregar<-sqlAppendTable(conexion_base, 'puntos_carga_local',guardar_datos)
        
        dbGetQuery(conexion_base,query_agregar)
        
        showNotification(
            h4("Datos guardados con éxito"), 
            action = NULL, duration = 5, type = "message")
    })
    
    ### Guarda los datos UTM
    observeEvent(input$inicio_correccion_UTM,{
        showModal(modalDialog(title="Guardar",
                              fluidRow(h4("Asgurate de haber elegido los nombres de las coordenadas correctamente antes de continuar")),
                              footer = tagList(
                                  actionButton("cancelar","Cancelar"),
                                  actionButton("aceptar_proceso_utm","Continuar")
                              )
        )) 
    })
    
    observeEvent(input$aceptar_proceso_utm,{
        # Restricciones
        if(length(input$tabla_inicio_utm_rows_selected)==0){
            showNotification(
                h4("Selecciona un renglón correspondiente con el punto pivote"), 
                action = NULL, duration = 5, type = "warning")
            
            removeModal()
            return()
        }
        guardar_input<-datos$puntos_coordenadasUTM[,c(input$col_nombre_UTM, input$col_N_UTM, input$col_E_UTM, input$col_h_UTM)]
        names(guardar_input)<-c("punto","n","e", "he")
        ### Restricciones
        if(class(guardar_input$n)!="numeric" || class(guardar_input$e)!="numeric" || 
           class(guardar_input$e)!="numeric" || class(guardar_input$punto)!="character"){
            showNotification(
                h4("Las columnas que seleccionaste no son numéricas; No se guardaran los datos"), 
                action = NULL, duration = 5, type = "warning")
            return()
        }
        
        user_id<-data.frame(user_id=rep(user_data()$user_id,length(datos$puntos_coordenadasUTM[,1])))
        fecha<-data.frame(fecha=rep(Sys.Date(),length(datos$puntos_coordenadasUTM[,1])))
        
        guardar_datos<-cbind(user_id, fecha, guardar_input)
        
        guardar_datos$fecha<-as.character(guardar_datos$fecha)
        query_agregar<-sqlAppendTable(conexion_base, 'puntos_utm',guardar_datos)
        
        dbGetQuery(conexion_base,query_agregar)
        
        
        showNotification(
            h4("Datos guardados con éxito"), 
            action = NULL, duration = 5, type = "message")
        
    })
    
    
})
