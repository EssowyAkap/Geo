# Dokumentacja
Przedstawiona dokumentacja opisuje działanie projektu zaliczeniowego nr 1 z przedmiotu ,,informatyka geodezyjna 2''. Celem programu jest transformowanie wprowadzanych (importowanych) współrzędnych z podanego układu do innego. Kalkulator geodezyjny został napisany przy pomocy języka programowania Python v3.11 i wymaga zainstalowania bibliotek numpy oraz argparse. Wykorzystane systemy operacyjne to Windows 10, a także macOS Sonoma 14.2.1 - działanie dla pozostałych systemów operacyjnych jeszcze nie zostało zweryfiowane.  
# Obsługiwane transformacje
XYZ (geocentryczne) -> BLH (elipsoidalne, fi, lambda)
BLH (elipsoidalne, fi, lambda) -> XYZ (geocentryczne)
XYZ (geocentryczne) -> NEUp (topocentryczne: north, east, up)
BL(GRS80, WGS84, ew. Krasowski) -> PL-2000 (w programie funkcja PL2000)
BL(GRS80, WGS84, ew. Krasowski) -> PL-1992 (w programie funkcja PL1992)

# Działanie programu
Program prosi o wywołanie takich rzeczy jak: nazwa elipsoidy, ścieżki do pliku tekstowego (który zawiera współrzędne, które użytkownik chce przetransformować), a także o rodzaj transformacji do wykonania. W przypadku odpowiednio wykonanych ,,inputów'' program wykonuje zadanie. Jeśli, któryś z kroków został wykonany nieprawidłowo, to użytkownik zostaje o tym w odpowiedni sposób poinformowany. 

# Instrukcja
Po uruchomieniu programu użytkownik zostanie poproszony o wybranie elipsoidy z pośród: "WGS84", "GRS80", "KRASOWSKI" - daje to kalkulatorowi informacje o parametrach danej elipsoidy. Następnie, użytkownik zostanie poproszony o podanie ścieżki do pliku .txt, który zawiera współrzędne do przekształcenia. Później należy wybrać rodzaj docelowej transformacji ("PL2000", "PL1992", "XYZ2NEUP", "XYZ2BLH", "BLH2XYZ"). Po zatwierdzeniu, program tworzy plik .txt (w lokalizacji, której on sam się znajduje), zawierający wyniki danych po przetransformowaniu. Użytkownik może wybrać, czy chce zakończyć działanie programu (klikając ENTER), czy dokonać kolejnego przeliczenia (równoznaczne z wpisaniem "DALEJ"). Program wymaga podania współrzędnych geocentrycznych w metrach, a elipsoidalnych w stopniach dziesiętnych. W takich samych jednostkach wyznacza te jednostki i zwraca użytkownikowi. Przykładowe użycie:

![image](https://github.com/EssowyAkap/Geo/assets/168012795/501dffbf-0436-4724-aac0-52aab42baa19)
1. Program wydaje polecenie ,,Elipsoida:''. Użytkownik musi wpisać nazwę danej elipsoidy, np. GRS80, po czym zatwierdza to przyciskiem ENTER.
![image](https://github.com/EssowyAkap/Geo/assets/168012795/7990f522-e9fa-4fa3-80d9-36232f13be9a)
2. Program prosi o ścieżkę do pliku z danymi, które ma wykorzstać w celu ich przeliczenia. Użystkownik podaje ścieżkę oraz potwierdza ją tym samym sposobem co w punkcie numer 1. Pojawia się prośba o podanie rodzaju transformacji:
![image](https://github.com/EssowyAkap/Geo/assets/168012795/2419a157-6ff1-4c6d-9d3e-f50c1a641cb0)
3. 



# Znane błędy i uwagi
Program zwraca informacje, o błędach popełnionych przez użytkownika - w przypadku podania niewłaściwego formatu danych w pliku .txt; braku przywoływanego pliku w danym folderze; błędnego wprowadzenia nazwy danego rodzaju transformaji lub elipsoidy. Dodatkowo, czasem przy użyciu transformacji "XYZ2NEUP" program tworzy pusty plik.
Dane w pliku .txt muszą być oddzielone od siebie spacją; dodatkowo muszą posiadać "." jako rozdzielenie części dziesiętnych od jedności.
Program pozwala na transformowanie współrzędnych ustawionych w kilku kolejnych wierszach. 

