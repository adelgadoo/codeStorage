// Initialize the library array
var myLibrary = [];

// Add a couple test books into the array
myLibrary[0] = "The Philosophy of WestWorld";
myLibrary[1] = "Metaphors and the words we live by";

// Set bookList as a global variable to append to
var bookList = document.getElementById("bookList");

function displayLibrary()
{ 
    // Loop through every existing book in library on load  
    for ( var i = 0; i < myLibrary.length; i++)
    {
        console.log(myLibrary[i]);
        var node = document.createElement("li");
        var bookToAdd = document.createTextNode(myLibrary[i]);
        node.appendChild(bookToAdd);
        bookList.appendChild(node);
    }
}

function addToLibrary()
{

    // Grab the book name from the text box
    var bookToAdd = document.getElementById("addBook").value;
    console.log(bookToAdd);

    // Capitalize the first word
    let bookToAddUpper = bookToAdd.replace(/^\w/, c=>c.toUpperCase());

    // Push the book name to the array of books in library
    myLibrary.push(bookToAddUpper);
    console.log(myLibrary);

    // Append to the list in the DOM
    addToList(bookToAddUpper);

    // Clear the text box
    document.getElementById("addBook").value = "";
}

function addToList(bookName)
{
    // Create the li element for the list
    let node = document.createElement("li");

    // Create text node for the book name
    let bookAddToList = document.createTextNode(bookName);

    // Append child
    node.appendChild(bookAddToList);
    bookList.appendChild(node);
}

function removeBook()
{
    
}

// Get the code to load 
window.onload = displayLibrary();

// Added code to listen for Enter to input a book
var input = document.getElementById("addBook");
input.addEventListener("keyup", function(event) {
    // Cancel the default action, if needed
    event.preventDefault();
    // Number 13 is the "Enter" key on the keyboard
    if (event.keyCode === 13) {
      // Trigger the button element with a click
      document.getElementById("addBtn").click();
    }
  });