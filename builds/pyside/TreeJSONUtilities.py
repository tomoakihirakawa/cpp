import json
import sys
from typing import Any, List, Dict, Union
from PySide6.QtWidgets import QTreeView, QApplication, QHeaderView
from PySide6.QtCore import QAbstractItemModel, QModelIndex, QObject, Qt, QFileInfo
from PySide6.QtGui import QBrush, QColor

class TreeItem:
    """A Json item corresponding to a line in QTreeView"""

    def __init__(self, parent: "TreeItem" = None):
        self._parent = parent
        self._key = ""
        self._value = ""
        self._value_type = None
        self._children = []

    def appendChild(self, item: "TreeItem"):
        """Add item as a child"""
        self._children.append(item)

    def child(self, row: int) -> "TreeItem":
        """Return the child of the current item from the given row"""
        return self._children[row]

    def parent(self) -> "TreeItem":
        """Return the parent of the current item"""
        return self._parent

    def childCount(self) -> int:
        """Return the number of children of the current item"""
        return len(self._children)

    def row(self) -> int:
        """Return the row where the current item occupies in the parent"""
        return self._parent._children.index(self) if self._parent else 0

    def toDict(self):
        """Recursively converts the TreeItem and its children into a dictionary or list."""
        if self.value_type is dict:
            return {child.key: child.toDict() for child in self._children}
        elif self.value_type is list:
            return [child.toDict() for child in self._children]
        else:
            return self.value

    @property
    def key(self) -> str:
        """Return the key name"""
        return self._key

    @key.setter
    def key(self, key: str):
        """Set key name of the current item"""
        self._key = key

    @property
    def value(self) -> str:
        """Return the value name of the current item"""
        return self._value

    @value.setter
    def value(self, value: str):
        """Set value name of the current item"""
        self._value = value

    @property
    def value_type(self):
        """Return the python type of the item's value."""
        return self._value_type

    @value_type.setter
    def value_type(self, value):
        """Set the python type of the item's value."""
        self._value_type = value

    @classmethod
    def load(cls, value: Union[List, Dict], parent: "TreeItem" = None, sort=True) -> "TreeItem":
        """Create a 'root' TreeItem from a nested list or a nested dictonary

        Examples:
            with open("file.json") as file:
                data = json.dump(file)
                root = TreeItem.load(data)

        This method is a recursive function that calls itself.

        Returns:
            TreeItem: TreeItem
        """
        rootItem = TreeItem(parent)
        rootItem.key = "root"

        if isinstance(value, dict):
            items = sorted(value.items()) if sort else value.items()

            for key, value in items:
                child = cls.load(value, rootItem)
                child.key = key
                child.value_type = type(value)
                rootItem.appendChild(child)

        elif isinstance(value, list):
            for index, value in enumerate(value):
                # delete filename from the child
                key = index
                if isinstance(value, dict) and "filename" in value:
                    key = value["filename"]
                    # del value["filename"]

                # ignore key of 'filename'
                child = cls.load(value, rootItem)
                # if there is filename, use it as key
                child.key = key


                child.value_type = type(value)
                rootItem.appendChild(child)

        else:
            rootItem.value = value
            rootItem.value_type = type(value)

        return rootItem

class JsonModel(QAbstractItemModel):
    """ An editable model of Json data """

    def __init__(self, parent: QObject = None):
        super().__init__(parent)

        self._rootItem = TreeItem()
        self._headers = ("key", "value")

    def toDict(self):
        """Converts the entire model's data into a dictionary or list based on the root item's type."""
        return self._rootItem.toDict()

    def clear(self):
        """ Clear data from the model """
        self.load({})

    def load(self, document: dict):
        """Load model from a nested dictionary returned by json.loads()

        Arguments:
            document (dict): JSON-compatible dictionary
        """

        assert isinstance(
            document, (dict, list, tuple)
        ), "`document` must be of dict, list or tuple, " f"not {type(document)}"

        self.beginResetModel()

        self._rootItem = TreeItem.load(document)
        self._rootItem.value_type = type(document)

        self.endResetModel()

        return True

    def data(self, index: QModelIndex, role: Qt.ItemDataRole) -> Any:
        if not index.isValid():
            return None

        item = index.internalPointer()

        # Display role for showing text
        if role == Qt.DisplayRole:
            if index.column() == 0:
                return item.key
            if index.column() == 1:
                return item.value

        # Edit role for editing value
        elif role == Qt.EditRole:
            if index.column() == 1:
                return item.value

        # Background role for changing background color based on type
        elif role == Qt.BackgroundRole:
            if hasattr(item, 'value_type') and item.value_type == dict and 'type' in item.value:
                type_value = item.value['type']
                if type_value == "RigidBody":
                    return QBrush(QColor(255, 229, 204))  # Light orange
                elif type_value == "Fluid":
                    return QBrush(QColor(204, 229, 255))  # Light blue
                # You can add more conditions for other types
            # Optionally, handle list items differently if needed

        return None

    def setData(self, index: QModelIndex, value: Any, role: Qt.ItemDataRole):
        """Override from QAbstractItemModel

        Set json item according index and role

        Args:
            index (QModelIndex)
            value (Any)
            role (Qt.ItemDataRole)

        """
        if role == Qt.EditRole:
            if index.column() == 1:
                item = index.internalPointer()
                item.value = str(value)

                self.dataChanged.emit(index, index, [Qt.EditRole])

                return True

        return False

    def headerData(
        self, section: int, orientation: Qt.Orientation, role: Qt.ItemDataRole
    ):
        """Override from QAbstractItemModel

        For the JsonModel, it returns only data for columns (orientation = Horizontal)

        """
        if role != Qt.DisplayRole:
            return None

        if orientation == Qt.Horizontal:
            return self._headers[section]

    def index(self, row: int, column: int, parent=QModelIndex()) -> QModelIndex:
        """Override from QAbstractItemModel

        Return index according row, column and parent

        """
        if not self.hasIndex(row, column, parent):
            return QModelIndex()

        if not parent.isValid():
            parentItem = self._rootItem
        else:
            parentItem = parent.internalPointer()

        childItem = parentItem.child(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QModelIndex()

    def parent(self, index: QModelIndex) -> QModelIndex:
        """Override from QAbstractItemModel

        Return parent index of index

        """

        if not index.isValid():
            return QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.parent()

        if parentItem == self._rootItem:
            return QModelIndex()

        return self.createIndex(parentItem.row(), 0, parentItem)

    def rowCount(self, parent=QModelIndex()):
        """Override from QAbstractItemModel

        Return row count from parent index
        """
        if parent.column() > 0:
            return 0

        if not parent.isValid():
            parentItem = self._rootItem
        else:
            parentItem = parent.internalPointer()

        return parentItem.childCount()

    def columnCount(self, parent=QModelIndex()):
        """Override from QAbstractItemModel

        Return column number. For the model, it always return 2 columns
        """
        return 2

    def flags(self, index: QModelIndex) -> Qt.ItemFlags:
        """Override from QAbstractItemModel

        Return flags of index
        """
        flags = super(JsonModel, self).flags(index)

        if index.column() == 1:
            return Qt.ItemIsEditable | flags
        else:
            return flags

    def to_json(self, item=None):

        if item is None:
            item = self._rootItem

        nchild = item.childCount()

        if item.value_type is dict:
            document = {}
            for i in range(nchild):
                ch = item.child(i)
                document[ch.key] = self.to_json(ch)
            return document

        elif item.value_type == list:
            document = []
            for i in range(nchild):
                ch = item.child(i)
                document.append(self.to_json(ch))
            return document

        else:
            return item.value

    def findItemByFilename(self, filename: str, parentItem: TreeItem = None) -> QModelIndex:
        """
        Recursively search for an item with the given filename (key) starting from parentItem.
        
        Args:
            filename (str): The filename (or key) to search for.
            parentItem (TreeItem): The item to start the search from. If None, start from the root.
        
        Returns:
            QModelIndex: The index of the found item, or an invalid QModelIndex if not found.
        """
        if parentItem is None:
            parentItem = self._rootItem
        
        # Check if the current item matches the filename
        if parentItem.key == filename:
            # We use 0 for the column as we're primarily interested in the item itself, not a specific column.
            return self.createIndex(parentItem.row(), 0, parentItem)
        
        # Recursively search in child items
        for i in range(parentItem.childCount()):
            childItem = parentItem.child(i)
            resultIndex = self.findItemByFilename(filename, childItem)
            if resultIndex.isValid():
                return resultIndex
        
        # If we reach here, the item was not found
        return QModelIndex()

    def findInInputFilesJSON(self, filename: str) -> QModelIndex:
        """
        Public method to start the search for an item with the given filename.
        
        Args:
            filename (str): The filename to search for.
        
        Returns:
            QModelIndex: The index of the found item, or an invalid index if not found.
        """
        return self.findItemByFilename(filename)


    # def data(self, index, role=Qt.DisplayRole):
    #     ItemTypeRole = Qt.UserRole + 1  # Custom role for item type
    #     if not index.isValid():
    #         return None
    #     item = index.internalPointer()
    #     if role == Qt.DisplayRole:
    #         if index.column() == 0:
    #             return item.key
    #         elif index.column() == 1:
    #             return item.value
    #     elif role == ItemTypeRole:
    #         # Assume each TreeItem has a 'type' attribute; adjust as necessary
    #         return getattr(item, 'type', None)
