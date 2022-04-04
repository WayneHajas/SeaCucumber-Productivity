import os
import sys
import win32com.client
import pythoncom

#This is a a revision to a previous class, DAO, that I created.
#DAO used Microsoft DAO 3.6 Object
#This class uses Microsoft ADO

#I expect backwards compatibility.  This class should work with files from ACCESS97 to ACCESS2010.

#Microsoft is no longer actively supportig DAO and therefore I am updating the class before there is a crisis.


#DAO taken from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/303349.
#Acknowledgement to Larry Bates
'''
20190826
   Removed references to pdb library
   import pythoncom and add CoInitialize command
      I suddenly started getting an error-message that said:  CoInitialize has not been called.
'''

class adoBaseClass:
    '''
    This base class is used to read/write data from/to daoDatabase.  It
    uses ADO    
    '''
    #
    # Completely free, no warranties.
    #
    def __init__(self, databasepath, SQL_query=None):
        '''adoBaseClass(databasepath, SQL_query=None)
        databasepath is the full name for the .mdb or .accdb file
        SQL_query is an optional query '''
        #----------------------------------------------------------------------
        #
        self.databasepath=databasepath        
        self.oRS=None #Will be given a value when an query is submitted
        self.Nrec=None
        #
        try:    
            pythoncom.CoInitialize()
            self.oConn = win32com.client.Dispatch('ADODB.Connection')
            self.oConn.ConnectionString = "Provider=Microsoft.ACE.OLEDB.12.0;Data source="+databasepath+";"
            self.oConn.Open()
            self.oRS = win32com.client.Dispatch('ADODB.Command')
            self.oRS.ActiveConnection = self.oConn    # Set the recordset to connect thru oConn
        except:
            print ('ADO line 40')
            print(("Unable to open ",databasepath))
            #retry the failed commands and hope for a helpful error message  
            self.oConn = win32com.client.Dispatch('ADODB.Connection')
            self.oConn.ConnectionString = "Provider=Microsoft.ACE.OLEDB.12.0;Data source="+databasepath+";Mode=16;"
            self.oConn.Open()
            self.oRS = win32com.client.Dispatch('ADODB.Command')
            self.oRS.ActiveConnection = self.oConn    # Set the recordset to connect thru oConn          
  
        #----------------------------------------------------------------------
        # If query is specified, open a recordset with it
        #----------------------------------------------------------------------
        if SQL_query!=None: self.execute(SQL_query)       
        return    
        
    def execute(self, SQL_query):
        '''
        ADO.execute(SQL_query)
        This method is used to execute a SQL_query against the    database.
        '''
        try:
            self.oRS = win32com.client.Dispatch('ADODB.Command')
            self.oRS.ActiveConnection = self.oConn    # Set the recordset to connect thru oConn
            self.oRS.CommandText = SQL_query        
        except:
            print('ADO line 65')
            print(('failed to execute\n',SQL_query))
            #repeat the failed operation incase there is a useful error message 
            self.oRS = win32com.client.Dispatch('ADODB.Command')
            self.oRS.ActiveConnection = self.oConn    # Set the recordset to connect thru oConn
            self.oRS.CommandText = SQL_query

        try:
            self.rs = self.oRS.Execute(1,1)[0]
        except:
            print('ADO line 75')
            print (SQL_query)
            print (self.oConn.ConnectionString)
            #repeat the failed operation incase there is a useful error message 
            self.rs =self.oRS.Execute(1,1)[0]

        self.DefineFieldNames()    
        try:
            self.GetNrec()
            self.rs.MoveFirst()
        except:
            self.Nrec=0
            

        return

    def __getitem__(self, key):
        #----------------------------------------------------------------------
        # Try to get the Value for the field requested from the recordset
        #----------------------------------------------------------------------
        try:    return self.rs.Fields(key).Value
        except:
            print('ado line 81')
            print(('Failure retreiving field named ',key))
            return self.rs.Fields(key).Value#In case there is a useful error message

    def __setitem__(self, key, value):
        '''ado._setitem(key,value)
        In the current record, set the value of the key-field to value'''
        
        try:
            self.rsFields.Item(key).Value=value
        except:
            print('ADO line 92')
            print(('unable to assign a value of ',value, ' to a field named ',key))
            self.rs.Fields.Item(key).Value=value#In case there is a helpful error message
            


    
    def MoveNext(self):
        '''ADO.MoveNext()
        Go to the next record'''
        self.rs.MoveNext()
        return



    def Update(self):
        '''ADO.Update
        Update values in current record'''
        self.rs.Update()
        return

    def close(self):
        self.rs.Close()


    def Fields(self):
        result=self.rs.Fields
        return result

    def GetFieldCount(self):
       result=len(self.rs.Fields)
       return result    
    def DefineFieldNames(self):
       'Get field names and put them in a list'
       try:
           self.Fname=[f.Name for f in self.rs.Fields]
       except:
           print('ADO line 134')
           self.Fname=[f.Name for f in self.rs.Fields] #in case there is a helpful error message
    
       return

    def Get(self):
        'Get Values from Current record and move to next record.  result in list'
        if (self.rs.EOF):return(None)
        result=[f.Value for f in self.rs.Fields]
        if not(self.rs.EOF ):
            self.rs.MoveNext()
        return(result)
        
    def GetNrec(self):
        self.Nrec=0
        if(self.rs.EOF):return
        try:
            self.rs.MoveFirst()
            self.Nrec=1
            while(not(self.rs.EOF )):
                self.Nrec+=1   
                self.rs.MoveNext()
        except:
            return
            
    
    def GetALL(self):
        'Get all values from the query'
        if self.Nrec==0:
            #nfield=self.rs.Fields
            #result=0*[nfield*None]
            #return(result)
            return([])
        self.rs.MoveFirst()
        result=[]
        while(not(self.rs.EOF )):
            result+=[self.Get()]
       
        return(result)
    
    def GetRec(self, RecNum):
        'Get a particular record'
        if (RecNum<0):return (None)
        if (RecNum>(.5+self.Nrec)):return (None)
        if (self.Nrec==0):return(None)
        self.rs.MoveFirst()
        for i in range(RecNum):
            result=self.Get()
        return(result)
    
    
    def GetVariable(self,name):
        if not name in self.Fname:
            return (None)
        if (self.Nrec==0):return([])
        self.rs.MoveFirst()
        result=[]
        while(not(self.rs.EOF )):
            result+=[self.rs.Fields(name).Value]
            self.rs.MoveNext()
        self.rs.MoveFirst()
        return(result)

            
if __name__ == "__main__":
    SQL_query='SELECT DISTINCT Headers.SurveyTitle as SurveyTitle, Headers.Year as yr FROM Headers WHERE ( ((Headers.SurveyTitle) Is Not Null)  AND ((Headers.Year) Is Not Null) ) ORDER BY Headers.SurveyTitle, Headers.Year;'
    SQL_query='SELECT  HeadersOri.SurveyTitle FROM HeadersOri;'
    databasepath='H:\AnalysisPrograms2013\PyFunctions\Geoduck\SampleAIP\Geoduck_BioNew.mdb'
    test1=adoBaseClass(databasepath)
    #print test1
    test2=adoBaseClass(databasepath,SQL_query)
    #test2.rs.MoveFirst()
    #print(test2.rs.EOF)
    test3=test2.GetALL()
    #print test
    #print test2.Nrec, ' Nrec'
    #print test3[test2.Nrec-1]
    print((test2.GetRec(1)))
    print((test2.GetVariable('Year')))
      
    
