# DrugSearch

---
This directory of the repository contains the code for the **DrugSearch** project.

## What is DrugSearch?

---
**TL;DR:** DrugSearch is a web application that is developed as a final year project at KCL.
It allows users to search for drugs and their interactions with other drugs.

**Full explanation:** DrugSearch is a web application that works with the **FDA** database which
contains information about drugs and their interactions with other drugs. The point of the application
is to visualize the data available in the database and make it easier for medical practitioners
to understand the adverse effects, interactions, and other information about the medications
they prescribe. Also, the application is targeted not only at physicians but also regular users
without any medical background.

## How is it developed?

---
The application is developed using **Flask**/**React** combo. The backend is handled by Flask
and the frontend is handled by React. Since the app will not be only sending and receiving data
through API calls bu also handling sensitive user data, this user data is going to be stored in a
secure way.
