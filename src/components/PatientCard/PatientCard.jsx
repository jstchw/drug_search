import React, {useEffect, useRef} from 'react';
import { Card, Badge } from 'react-bootstrap';
import { searchSex } from "../OptionModal/OptionModal";


const processSearchOptions = (searchOptions, setSearchOptions) => {

    let newSearchOptions = {...searchOptions}

    if(searchOptions.age[0] && searchOptions.age[1]) {
        console.log(searchOptions.age)
        newSearchOptions = {
            ...newSearchOptions,
            age: searchOptions.age.join(' - ')
        }
    }
    // If value from searchSex is equal to the value of searchOptions.sex then return the label
    if(searchOptions.sex) {
        newSearchOptions = {
            ...newSearchOptions,
            sex: searchSex.find(sex => sex.value === searchOptions.sex).label
        }
    }

    setSearchOptions(newSearchOptions)
}

const PatientCard = (props) => {
    const { searchBy, ...rest } = props.searchOptions;
    const [searchOptions, setSearchOptions] = React.useState(rest)

    useEffect(() => {
        processSearchOptions(searchOptions, setSearchOptions)
    }, []);

    return (
        <Card className={'my-3'}>
            <Card.Header>
                <span>Patient Card</span>
            </Card.Header>
            <Card.Body className={'text-md-start'}>
                {
                    Object.keys(searchOptions).map((key, index) => {
                        return (
                            <Card.Text key={index}>
                                <Badge>
                                    {
                                        key.charAt(0).toUpperCase() + key.slice(1)
                                    }
                                </Badge>
                                &nbsp;
                                <span>
                                    {
                                        searchOptions[key]
                                    }
                                </span>
                            </Card.Text>
                        )
                    })
                }
            </Card.Body>
        </Card>
    )
}

export default PatientCard