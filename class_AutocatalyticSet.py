"""
Author: Dr. techn. Sebastian Raubitzek MSc. BSc.; SBA Research, Complexity and Resilience Research Group
The AutocatalyticSet class models and analyzes autocatalytic reaction networks, focusing on the synthesis of complex products from simpler constituents. This class is central to the study presented in "Autocatalytic Sets and Assembly Theory: A Toy Model Perspective" by Sebastian Raubitzek et al., where the emergence of complex structures through autocatalysis is explored.

### Purpose:
The AutocatalyticSet class simulates the formation of a target product (final product) through a series of catalytic reactions, where the products of some reactions act as catalysts for others. The class allows for the exploration of different catalyst assignment strategies and provides tools for visualizing and analyzing the resulting reaction networks.

### Key Features:
- **Initialization**: The class is initialized with a final product, catalyst assignment probability, and a strategy for how catalysts are assigned.
- **Catalyst Assignment Strategies**:
  - **Random assignment without reuse** (`switch_cat = 0`): Catalysts are assigned randomly to reactions, with no reuse of catalysts.
  - **Weighted random assignment without reuse** (`switch_cat = 1`): Catalysts are assigned based on a weighted probability proportional to their complexity, without reuse.
  - **Weighted random assignment with reuse** (`switch_cat = 2`): Catalysts are assigned based on weighted probabilities, with the possibility of reuse across multiple reactions.
- **Reaction Network Construction**: The class recursively splits the final product into smaller reactants, simulating a hierarchical assembly process and generating a network of reactions.
- **Visualization**: The reaction network can be visualized using NetworkX and Matplotlib, with customizable color mappings for different elements.
- **Subnetwork Analysis**: The class can extract and analyze subnetworks within the larger reaction network, checking for autocatalytic properties and allowing detailed study of smaller reaction clusters.

### Attributes:
- **final_product (str)**: The final product to be synthesized through the reactions.
- **catalyst_probability (float)**: The probability of assigning a catalyst to each reaction, ranging from 0 to 1.
- **switch_cat (int)**: Determines the strategy used for catalyst assignment.
- **reactions (list)**: A list that stores the reactions as they are generated.
- **reactions_df (pd.DataFrame)**: A pandas DataFrame that holds the reaction data.
- **elements (list)**: A list of elements (constituents) in the final product.
- **color_map (dict)**: A dictionary mapping elements to colors for visualization purposes.
- **subnetwork_dict (dict)**: A dictionary storing extracted subnetworks for further analysis.

### Methods:
- **assign_catalysts()**: Assigns catalysts to each reaction based on the selected strategy and probability.
- **calculate_weighted_probabilities(elements)**: Calculates the weighted probabilities for elements based on their complexity (length).
- **assign_colors(elements)**: Assigns unique colors to elements for visualization.
- **add_reaction(layer, reactant1, reactant2, product, catalyst)**: Adds a reaction to the network, combining two reactants into a product.
- **split_product(product, layer)**: Recursively splits the final product into smaller reactants, generating the reaction network.
- **build_reactions()**: Constructs the entire reaction network starting from the final product.
- **save_reactions_to_file(filename)**: Saves the reactions to a CSV file for external analysis.
- **display_reactions()**: Displays the reactions in a tabular format.
- **draw_network()**: Visualizes the entire reaction network.
- **get_layers(G)**: Determines the hierarchical layers of the reaction network for visualization.
- **draw_network_sub(reactions_df, key)**: Visualizes a subnetwork of the reaction network.
- **get_layers_sub(G, key)**: Determines the hierarchical layers of a subnetwork for visualization.
- **extract_subnetworks()**: Extracts subnetworks from the main reaction network for detailed analysis.
- **is_autocatalytic_subnetwork(subnetwork_df, check_for_catalysts=True)**: Checks if a given subnetwork is autocatalytic.

### Context and Usage:
The AutocatalyticSet class is designed to be used within the broader framework of the study on autocatalytic sets and assembly theory. It enables the simulation and analysis of reaction networks, providing insights into how complex systems can self-organize through autocatalysis. The class can be used in conjunction with other tools and scripts provided in the project to perform extensive analyses, including comparing different reaction networks, analyzing their properties, and visualizing the assembly pathways.

This class is particularly valuable for exploring theoretical concepts related to autocatalysis and understanding how small changes in catalyst assignment or reaction probabilities can lead to significant differences in the assembly of complex structures. The visualizations generated by the class help in intuitively grasping the underlying processes and dynamics of these networks.
"""

import numpy as np
import random
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class AutocatalyticSet:
    def __init__(self, final_product, catalyst_probability=1.0, switch_cat=0, print_out=True):
        """
        Initializes the AutocatalyticSet with the final product and catalyst assignment strategy.

        Args:
        final_product (str): The final product to be synthesized through the reactions.
        catalyst_probability (float): Probability of assigning a catalyst to each reaction (0 to 1).
        switch_cat (int): Switch for catalyst assignment strategy (0: random, 1: weighted random, 2: weighted random with reuse).
        """
        self.print_out = print_out
        self.final_product = final_product
        self.reactions = []  # List to store reactions
        self.reactions_df = None
        self.elements = list(final_product)  # List of elements in the final product
        self.catalyst_probability = catalyst_probability  # Probability of assigning a catalyst
        self.switch_cat = switch_cat  # Catalyst assignment strategy
        self.color_map = self.assign_colors(self.elements)  # Map elements to colors for visualization
        self.subnetwork_dict = {}  # Dictionary to store subnetworks

    def assign_catalysts(self):
        """
        Assigns catalysts to each reaction based on the defined probability and strategy.
        """
        all_elements = set(self.elements)
        all_products = {reaction['Product'] for reaction in self.reactions}
        all_possible_catalysts = list(all_elements | all_products)
        used_catalysts = set()

        if self.switch_cat == 0:
            # Random assignment without reuse
            for reaction in self.reactions:
                if random.random() <= self.catalyst_probability:
                    possible_catalysts = [e for e in all_possible_catalysts if
                                          e not in [reaction['Constituent 1'], reaction['Constituent 2']] and e not in used_catalysts]
                    if possible_catalysts:
                        reaction['Catalyst'] = random.choice(possible_catalysts)
                        used_catalysts.add(reaction['Catalyst'])
                    else:
                        reaction['Catalyst'] = ''
                else:
                    reaction['Catalyst'] = ''
        elif self.switch_cat == 1:
            # Weighted random assignment without reuse
            weighted_catalysts = self.calculate_weighted_probabilities(all_possible_catalysts)
            for reaction in self.reactions:
                if random.random() <= self.catalyst_probability:
                    possible_catalysts = [e for e in all_possible_catalysts if
                                          e not in [reaction['Constituent 1'], reaction['Constituent 2']] and e not in used_catalysts]
                    if possible_catalysts:
                        chosen_catalyst = random.choices(possible_catalysts, weights=[weighted_catalysts[c] for c in possible_catalysts])[0]
                        reaction['Catalyst'] = chosen_catalyst
                        used_catalysts.add(reaction['Catalyst'])
                    else:
                        reaction['Catalyst'] = ''
                else:
                    reaction['Catalyst'] = ''
        elif self.switch_cat == 2:
            # Weighted random assignment with reuse
            weighted_catalysts = self.calculate_weighted_probabilities(all_possible_catalysts)
            for reaction in self.reactions:
                if random.random() <= self.catalyst_probability:
                    possible_catalysts = [e for e in all_possible_catalysts if
                                          e not in [reaction['Constituent 1'], reaction['Constituent 2']]]
                    if possible_catalysts:
                        chosen_catalyst = random.choices(possible_catalysts, weights=[weighted_catalysts[c] for c in possible_catalysts])[0]
                        reaction['Catalyst'] = chosen_catalyst
                    else:
                        reaction['Catalyst'] = ''
                else:
                    reaction['Catalyst'] = ''

    def calculate_weighted_probabilities(self, elements):
        """
        Calculate weighted probabilities based on the length of the elements.
        The probability of an element being chosen as a catalyst is proportional to its length.

        Args:
        elements (list): List of elements to calculate probabilities for.

        Returns:
        dict: Mapping of elements to their calculated probabilities.
        """
        max_length = max(len(e) for e in elements)
        probabilities = {e: len(e) / max_length for e in elements}
        total_weight = sum(probabilities.values())
        normalized_probabilities = {e: p / total_weight for e, p in probabilities.items()}
        return normalized_probabilities

    def assign_colors(self, elements):
        """
        Assigns unique colors to each element for visualization purposes.

        Args:
        elements (list): List of elements to assign colors to.

        Returns:
        dict: Mapping of elements to their assigned colors.
        """
        colors = [
            'aliceblue', 'antiquewhite', 'black', 'blue', 'blueviolet', 'brown', 'cadetblue', 'chartreuse', 'chocolate',
            'coral', 'cornflowerblue', 'crimson', 'cyan', 'darkblue', 'darkcyan', 'darkgoldenrod', 'darkgray',
            'darkgreen', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid', 'darkred', 'darksalmon',
            'darkseagreen', 'darkslateblue', 'darkslategray', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue',
            'dimgray', 'dodgerblue', 'firebrick', 'forestgreen', 'fuchsia', 'gold', 'goldenrod', 'gray', 'green',
            'greenyellow', 'indianred', 'indigo', 'lawngreen', 'lightcyan', 'lightseagreen', 'lightslategray', 'lime',
            'limegreen', 'magenta', 'maroon', 'mediumblue', 'mediumorchid', 'mediumpurple', 'mediumseagreen',
            'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'navy', 'olive',
            'olivedrab', 'orange', 'orangered', 'orchid', 'palevioletred', 'papayawhip', 'peru', 'plum', 'purple',
            'rebeccapurple', 'red', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown', 'seagreen', 'sienna',
            'silver', 'skyblue', 'slateblue', 'slategray', 'springgreen', 'steelblue', 'tan', 'teal', 'thistle', 'tomato',
            'turquoise', 'violet', 'yellowgreen'
        ]
        random.shuffle(colors)
        return {element: colors[i % len(colors)] for i, element in enumerate(elements)}

    def add_reaction(self, layer, reactant1, reactant2, product, catalyst):
        """
        Adds a reaction to the set, combining two reactants into a product with a catalyst.

        Args:
        layer (int): The hierarchical layer of the reaction.
        reactant1 (str): The first reactant.
        reactant2 (str): The second reactant.
        product (str): The product of the reaction.
        catalyst (str): The catalyst for the reaction.
        """
        reaction_number = len(self.reactions)
        self.reactions.append({
            'Reaction': reaction_number,
            'Layer': layer,
            'Constituent 1': reactant1,
            'Constituent 2': reactant2,
            'Catalyst': catalyst,
            'Product': product
        })
        for element in [reactant1, reactant2, product, catalyst]:
            if element not in self.color_map:
                self.color_map[element] = random.choice(list(mcolors.CSS4_COLORS.keys()))

    def split_product(self, product, layer):
        """
        Recursively splits the product into smaller reactants and records the reactions.

        Args:
        product (str): The product to be split.
        layer (int): The hierarchical layer of the reaction.
        """
        if len(product) == 1:
            return
        split_index = random.randint(1, len(product) - 1)
        reactant1 = ''.join(sorted(product[:split_index]))
        reactant2 = ''.join(sorted(product[split_index:]))
        self.split_product(reactant1, layer + 1)
        self.split_product(reactant2, layer + 1)
        self.add_reaction(layer, reactant1, reactant2, product, catalyst='')

    def build_reactions(self):
        """
        Builds the full set of reactions starting from the final product.
        """
        self.split_product(self.final_product, 0)
        self.assign_catalysts()
        self.reactions_df = pd.DataFrame(self.reactions)

    def save_reactions_to_file(self, filename):
        """
        Saves the reactions to a CSV file.

        Args:
        filename (str): The name of the file to save the reactions.
        """
        df = pd.DataFrame(self.reactions)
        df.to_csv(filename, index=False)

    def display_reactions(self):
        """
        Prints the reactions as a DataFrame.
        """
        df = pd.DataFrame(self.reactions)
        print(df)

    def draw_network(self):
        """
        Visualizes the reaction network using NetworkX and Matplotlib.
        """
        G = nx.DiGraph()
        pos = {}
        labels = {}

        for reaction in self.reactions:
            G.add_edge(reaction['Constituent 1'], reaction['Product'], catalyst=reaction['Catalyst'])
            G.add_edge(reaction['Constituent 2'], reaction['Product'], catalyst=reaction['Catalyst'])
            labels[(reaction['Constituent 1'], reaction['Product'])] = reaction['Catalyst']
            labels[(reaction['Constituent 2'], reaction['Product'])] = reaction['Catalyst']

        layers = self.get_layers(G)
        num_layers = len(layers)

        for layer, nodes in enumerate(layers):
            sorted_nodes = sorted(nodes)
            for i, node in enumerate(sorted_nodes):
                pos[node] = (i - len(nodes) / 2.0, num_layers - layer - 1)

        node_colors = [self.color_map[node] for node in G.nodes()]
        plt.figure(figsize=(12, 8))
        nx.draw(G, pos, with_labels=True, node_size=3000, node_color=node_colors, font_size=10, arrowsize=20)

        edge_labels = nx.get_edge_attributes(G, 'catalyst')
        for edge, catalyst in edge_labels.items():
            nx.draw_networkx_edge_labels(G, pos, edge_labels={edge: catalyst}, font_color=self.color_map.get(catalyst, 'black'))

        for edge in G.edges(data=True):
            catalyst = edge[2]['catalyst']
            edge_color = self.color_map.get(catalyst, 'black') if catalyst else 'black'
            nx.draw_networkx_edges(G, pos, edgelist=[(edge[0], edge[1])], edge_color=edge_color, width=2.0)

        nx.draw_networkx_labels(G, pos, labels={node: node for node in G.nodes()},
                                font_color='white',
                                bbox=dict(facecolor='black', edgecolor='none', boxstyle='round,pad=0.5'))

        plt.savefig('reaction_network.png', format='png')
        plt.savefig('reaction_network.svg', format='svg')
        plt.show()

    def get_layers(self, G):
        """
        Determines the hierarchical layers of the network for visualization.

        Args:
        G (networkx.DiGraph): The directed graph representing the reaction network.

        Returns:
        list: A list of sets, each containing the nodes at a particular layer.
        """
        layers = []
        current_layer = {self.final_product}
        while current_layer:
            layers.append(current_layer)
            next_layer = set()
            for node in current_layer:
                predecessors = list(G.predecessors(node))
                next_layer.update(predecessors)
            current_layer = next_layer
        return layers

    def draw_network_sub(self, reactions_df, key):
        """
        Visualizes a subnetwork of the reaction network using NetworkX and Matplotlib.

        Args:
        reactions_df (pd.DataFrame): The dataframe containing reactions.
        key (str): The key node to visualize the subnetwork from.
        """
        G = nx.DiGraph()
        pos = {}
        labels = {}

        for _, reaction in reactions_df.iterrows():
            G.add_edge(reaction['Constituent 1'], reaction['Product'], catalyst=reaction['Catalyst'])
            G.add_edge(reaction['Constituent 2'], reaction['Product'], catalyst=reaction['Catalyst'])
            labels[(reaction['Constituent 1'], reaction['Product'])] = reaction['Catalyst']
            labels[(reaction['Constituent 2'], reaction['Product'])] = reaction['Catalyst']

        layers = self.get_layers_sub(G, key)
        num_layers = len(layers)

        for layer, nodes in enumerate(layers):
            sorted_nodes = sorted(nodes)
            for i, node in enumerate(sorted_nodes):
                pos[node] = (i - len(nodes) / 2.0, num_layers - layer - 1)

        node_colors = [self.color_map[node] for node in G.nodes()]
        plt.figure(figsize=(12, 8))
        nx.draw(G, pos, with_labels=True, node_size=3000, node_color=node_colors, font_size=10, arrowsize=20)

        edge_labels = nx.get_edge_attributes(G, 'catalyst')
        for edge, catalyst in edge_labels.items():
            nx.draw_networkx_edge_labels(G, pos, edge_labels={edge: catalyst}, font_color=self.color_map.get(catalyst, 'black'))

        for edge in G.edges(data=True):
            catalyst = edge[2]['catalyst']
            edge_color = self.color_map.get(catalyst, 'black') if catalyst else 'black'
            nx.draw_networkx_edges(G, pos, edgelist=[(edge[0], edge[1])], edge_color=edge_color, width=2.0)

        nx.draw_networkx_labels(G, pos, labels={node: node for node in G.nodes()},
                                font_color='white',
                                bbox=dict(facecolor='black', edgecolor='none', boxstyle='round,pad=0.2'))

        plt.show()

    def get_layers_sub(self, G, key):
        """
        Determines the hierarchical layers of a subnetwork for visualization.

        Args:
        G (networkx.DiGraph): The directed graph representing the reaction network.
        key (str): The key node to determine the layers from.

        Returns:
        list: A list of sets, each containing the nodes at a particular layer.
        """
        layers = []
        current_layer = {key}
        while current_layer:
            layers.append(current_layer)
            next_layer = set()
            for node in current_layer:
                predecessors = list(G.predecessors(node))
                next_layer.update(predecessors)
            current_layer = next_layer
        return layers

    def extract_subnetworks(self):
        """
        Extracts subnetworks from the main reaction dataframe and stores them in new dataframes.
        """

        def get_subnetwork(product):
            subnetwork = self.reactions_df[self.reactions_df['Product'] == product]
            reactants = set(subnetwork['Constituent 1'].tolist() + subnetwork['Constituent 2'].tolist())

            while reactants:
                next_reactants = set()
                for reactant in reactants:
                    sub_reactions = self.reactions_df[self.reactions_df['Product'] == reactant]
                    subnetwork = pd.concat([subnetwork, sub_reactions])
                    next_reactants.update(sub_reactions['Constituent 1'].tolist() + sub_reactions['Constituent 2'].tolist())
                reactants = next_reactants - set(subnetwork['Product'].tolist())

            # Remove networks with product == self.final_product, i.e. the whole network must not be part of the subnetworks
            #subnetwork = subnetwork[subnetwork['Product'] != self.final_product]

            return subnetwork.drop_duplicates().reset_index(drop=True)

        def process_product_print_out(product):
            subnetwork_df = get_subnetwork(product)
            if not subnetwork_df.empty:
                if product != self.final_product:
                    self.subnetwork_dict[product] = subnetwork_df
                    print(f"Subnetwork for {product}:")
                    print(subnetwork_df)
                # self.subnetwork_dict[product] = subnetwork_df
                constituents = set(subnetwork_df['Constituent 1'].tolist() + subnetwork_df['Constituent 2'].tolist())
                for constituent in constituents:
                    if constituent not in self.subnetwork_dict:
                        process_product_print_out(constituent)

        def process_product(product):
            subnetwork_df = get_subnetwork(product)
            if not subnetwork_df.empty:
                if product != self.final_product:
                    self.subnetwork_dict[product] = subnetwork_df
                #self.subnetwork_dict[product] = subnetwork_df
                constituents = set(subnetwork_df['Constituent 1'].tolist() + subnetwork_df['Constituent 2'].tolist())
                for constituent in constituents:
                    if constituent not in self.subnetwork_dict:
                        process_product(constituent)

        if self.print_out:
            process_product_print_out(self.final_product)
        else:
            process_product(self.final_product)


    def is_autocatalytic_subnetwork(self, subnetwork_df, check_for_catalysts=True):
        """
        Checks if a given subnetwork is autocatalytic.

        Args:
        subnetwork_df (pd.DataFrame): The dataframe representing the subnetwork.
        check_for_catalysts (bool): Whether to check if the subnetwork has at least one catalyst.

        Returns:
        bool: True if the subnetwork is autocatalytic, False otherwise.
        """
        products = set(subnetwork_df['Product'])
        constituents_1 = set(subnetwork_df['Constituent 1'])
        constituents_2 = set(subnetwork_df['Constituent 2'])
        all_elements = products | constituents_1 | constituents_2
        catalysts = set(subnetwork_df['Catalyst'].dropna())
        catalysts.discard('')
        if check_for_catalysts and not catalysts:
            return False
        if catalysts - all_elements:
            return False
        return True