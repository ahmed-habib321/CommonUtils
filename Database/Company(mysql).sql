CREATE DATABASE IF NOT EXISTS Company DEFAULT CHARACTER SET UTF8 DEFAULT COLLATE UTF8_GENERAL_CI;
USE Company;

CREATE TABLE IF NOT EXISTS Department (
        DepartmentNumber INT PRIMARY KEY,
        DepartmentName VARCHAR(150) UNIQUE,
        Location VARCHAR(100)
);

CREATE TABLE IF NOT EXISTS Project (
        ProjectNumber INT PRIMARY KEY,
        ProjectName VARCHAR(150) UNIQUE,
        Location VARCHAR(100),
        DepartmentNumber INT,
        FOREIGN KEY (DepartmentNumber) REFERENCES Department(DepartmentNumber)
);

CREATE TABLE IF NOT EXISTS Employee (
        EmployeeNumber INT PRIMARY KEY,
        EmployeeName VARCHAR(150) UNIQUE,
        Address VARCHAR(200),
        Salary DOUBLE,
        HiringDate DATE,
        BirthDate DATE,
        DepartmentNumber INT,
        FOREIGN KEY (DepartmentNumber) REFERENCES Department(DepartmentNumber)
);

CREATE TABLE IF NOT EXISTS EmployeePhones (
        Phone VARCHAR(20),
        EmployeeNumber INT,
        PRIMARY KEY (Phone, EmployeeNumber),
        FOREIGN KEY (EmployeeNumber) REFERENCES Employee(EmployeeNumber)
);

CREATE TABLE IF NOT EXISTS WorkOn (
        EmployeeNumber INT,
        ProjectNumber INT,
        PRIMARY KEY (ProjectNumber, EmployeeNumber),
        FOREIGN KEY (ProjectNumber) REFERENCES Project(ProjectNumber),
        FOREIGN KEY (EmployeeNumber) REFERENCES Employee(EmployeeNumber)
);

CREATE TABLE IF NOT EXISTS Users (
    Username VARCHAR(150) PRIMARY KEY,
    Password VARCHAR(150),
    IsBanned varchar(150)
);

CREATE VIEW IF NOT EXISTS EmployeeData AS 
        SELECT 
                EmployeeNumber AS 'Number',
                EmployeeName AS 'Name',
                Address,
                Salary,
                HiringDate AS 'HiringDate',
                BirthDate AS 'BirthDate',
                Employee.DepartmentNumber,
                DepartmentName AS 'Department'
        FROM 
                Department, Employee
        WHERE 
                Employee.DepartmentNumber = Department.DepartmentNumber;

CREATE VIEW IF NOT EXISTS DepartmentData AS
        SELECT 
                 DepartmentNumber,
                 DepartmentName AS 'Department',
                 Location
        FROM 
                Department;


CREATE VIEW IF NOT EXISTS ProjectData AS
        SELECT 
                ProjectNumber,
                ProjectName AS 'Project',
                Project.DepartmentNumber AS 'DepartmentNumber',
                Department.DepartmentName AS 'Department',
                Project.Location
        FROM 
                Project, Department
        WHERE 
                Project.DepartmentNumber = Department.DepartmentNumber;

CREATE VIEW IF NOT EXISTS EmployeePhonesData AS
        SELECT 
                Employee.EmployeeNumber AS 'EmployeeNumber',
                Employee.EmployeeName AS 'EmployeeName',
                EmployeePhones.Phone
        FROM 
                Employee, EmployeePhones
        WHERE 
                Employee.EmployeeNumber = EmployeePhones.EmployeeNumber; 

CREATE VIEW IF NOT EXISTS WorkOnData AS 
        SELECT 
                WorkOn.EmployeeNumber AS 'EmployeeNumber',
                Employee.EmployeeName AS 'EmployeeName',
                WorkOn.ProjectNumber AS 'ProjectNumber',
                Project.ProjectName AS 'ProjectName'
        FROM 
                WorkOn, Employee, Project
        WHERE 
                WorkOn.EmployeeNumber = Employee.EmployeeNumber 
                AND 
                WorkOn.ProjectNumber = Project.ProjectNumber;